import os
import subprocess
import tempfile
import textwrap
from pathlib import Path
from dotenv import load_dotenv

# 1. Environment and Path Setup
load_dotenv()
INPUT_DIR = os.getenv("INPUT_DIR", "../input")
REFERENCE_DIR = Path(INPUT_DIR) / "Reference"
REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

def zenodo_download(record_id, target_dir):
    """Download a Zenodo record using the zenodo-get CLI tool."""
    print(f"--- Downloading Zenodo record {record_id} to {target_dir} ---")
    try:
        # zenodo_get downloads to the current working directory
        subprocess.run(["zenodo_get", str(record_id)], cwd=str(target_dir), check=True)
        print(f"Download complete.")
        return True
    except FileNotFoundError:
        print("\nError: 'zenodo_get' command not found. Please install it: pip install zenodo-get")
        return False
    except subprocess.CalledProcessError as e:
        print(f"\nError downloading Zenodo record {record_id}: {e}")
        return False

def convert_rds_to_h5ad():
    """
    Uses SeuratDisk in R to convert the downloaded RDS file to h5ad.
    """
    print("\n--- Converting RDS to h5ad ---")
    
    # 1. Find the RDS file (zenodo_get might download scAtlas.rds.gz)
    potential_files = list(REFERENCE_DIR.glob("*.rds*"))
    if not potential_files:
        print(f"Error: No .rds files found in {REFERENCE_DIR}")
        return False
    
    rds_path = potential_files[0]
    h5seurat_path = rds_path.with_suffix('.h5seurat')
    h5ad_path = rds_path.with_suffix('.h5ad')
    
    # Standardize final name for the pipeline
    final_h5ad = REFERENCE_DIR / "loveless_2025_atlas.h5ad"

    print(f"Input file: {rds_path.name}")
    
    # 2. Generate the R script
    r_script_content = textwrap.dedent(f"""
        if (!requireNamespace("SeuratDisk", quietly = TRUE)) {{
            if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos='http://cran.us.r-project.org')
            remotes::install_github("mojaveazure/seurat-disk")
        }}
        
        library(SeuratDisk)
        
        print("Loading RDS file...")
        atlas <- readRDS("{rds_path.absolute()}")
        
        print("Saving as h5seurat...")
        SaveH5Seurat(atlas, filename = "{h5seurat_path.absolute()}", overwrite = TRUE)
        
        print("Converting to h5ad...")
        Convert("{h5seurat_path.absolute()}", dest = "h5ad", overwrite = TRUE)
        
        # Cleanup intermediate
        if (file.exists("{h5seurat_path.absolute()}")) {{
            file.remove("{h5seurat_path.absolute()}")
        }}
        
        print("Conversion complete!")
    """)

    with tempfile.NamedTemporaryFile(suffix=".R", mode="w", delete=False) as f:
        f.write(r_script_content)
        temp_r_script = f.name

    try:
        print("Executing R conversion (SeuratDisk)...")
        process = subprocess.Popen(["Rscript", temp_r_script], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end='')
        process.wait()
        
        if process.returncode == 0:
            # Rename to the expected pipeline name
            generated_h5ad = REFERENCE_DIR / (rds_path.stem + ".h5ad")
            if generated_h5ad.exists():
                generated_h5ad.rename(final_h5ad)
                print(f"Success: {final_h5ad}")
            elif h5ad_path.exists():
                h5ad_path.rename(final_h5ad)
                print(f"Success: {final_h5ad}")
            return True
        else:
            print(f"R script failed with return code {process.returncode}")
            return False
    except Exception as e:
        print(f"An error occurred: {e}")
        return False
    finally:
        if os.path.exists(temp_r_script):
            os.remove(temp_r_script)

def main():
    record_id = "14199536"
    
    if zenodo_download(record_id, REFERENCE_DIR):
        convert_rds_to_h5ad()

if __name__ == "__main__":
    main()
