import scanpy as sc
import squidpy as sq
import pandas as pd
from pathlib import Path
from MetScore import RankAndScore

def analyze_sample(sample_id, samples_dir="/scratch/rprest2/Spatial-MetScore/input/Samples", metadata_file="/scratch/rprest2/Spatial-MetScore/input/patient_metadata.csv", 
                   high_genes_file="/scratch/rprest2/Spatial-MetScore/input/met-high-genes-human.csv", low_genes_file="/scratch/rprest2/Spatial-MetScore/input/met-low-genes-human.csv",
                   min_counts=100, target_sum=1e4, load_images=False):
    """
    Process a single Visium spatial transcriptomics sample.
    
    Loads the sample data, filters and normalizes it, calculates the MetScore, 
    and retrieves the corresponding patient metadata.
    
    Parameters
    ----------
    sample_id : str
        The ID of the sample to process (e.g., "PT-1A"). Must match the folder name in `samples_dir`.
    samples_dir : str or Path
        Directory containing the sample folders.
    metadata_file : str or Path
        Path to the CSV file containing patient metadata.
    high_genes_file : str or Path
        Path to the CSV file containing the Met-high gene signature.
    low_genes_file : str or Path
        Path to the CSV file containing the Met-low gene signature.
    min_counts : int, optional
        Minimum number of counts required for a spot to be kept (default is 100).
    target_sum : float, optional
        Target sum for total count normalization (default is 1e4).
        
    Returns
    -------
    tuple
        A tuple containing:
        - AnnData: The processed AnnData object with calculated MetScores.
        - dict: A dictionary containing the sample's metadata.
    """
    
    sample_path = Path(samples_dir) / sample_id
    
    if not sample_path.exists():
        raise FileNotFoundError(f"Sample directory not found: {sample_path}")
        
    print(f"--- Processing Sample: {sample_id} ---")
    
    # 1. Load Metadata
    metadata = {}
    try:
        md_df = pd.read_csv(metadata_file)
        # Search for the row matching this Sample_ID (case-insensitive)
        match = md_df[md_df['Sample_ID'].str.upper() == sample_id.upper()]
        if not match.empty:
            metadata = match.iloc[0].to_dict()
            print(f"Loaded metadata for {sample_id}: {metadata.get('Tissue_Type', 'Unknown Tissue')}")
        else:
            print(f"Warning: No metadata found for {sample_id} in {metadata_file}")
    except Exception as e:
        print(f"Warning: Could not load metadata: {e}")

    # 2. Load Gene Signatures
    try:
        met_high_genes = pd.read_csv(high_genes_file)
        met_low_genes = pd.read_csv(low_genes_file)
        
        met_high_symbols = met_high_genes['SYMBOL'].tolist()
        met_low_symbols = met_low_genes['SYMBOL'].tolist()
    except Exception as e:
        raise RuntimeError(f"Failed to load gene signatures: {e}")

    # 3. Load Spatial Data
    print(f"Loading spatial data from {sample_path}...")
    try:
        adata = sq.read.visium(
            path=str(sample_path), 
            counts_file="filtered_feature_bc_matrix.h5",
            load_images=load_images
        )
        adata.var_names_make_unique()
        print(f"Loaded AnnData: {adata.n_obs} spots x {adata.n_vars} genes")
    except Exception as e:
        raise RuntimeError(f"Failed to load spatial data for {sample_id}: {e}")

    # 4. Preprocessing
    print("Preprocessing data (filtering, normalizing, log-transforming)...")
    sc.pp.filter_cells(adata, min_counts=min_counts)
    
    adata.layers['counts'] = adata.X.copy()
    
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    
    adata.layers['lognorm'] = adata.X.copy()

    # 5. Calculate MetScore
    print("Calculating MetScore...")
    adata = RankAndScore(adata, met_high_symbols, met_low_symbols)
    
    print("--- Processing Complete ---")
    return adata, metadata

if __name__ == "__main__":
    # Example usage:
    # adata, metadata = analyze_sample("PT-1A")
    # sq.pl.spatial_scatter(adata, color='MetScore', cmap='RdBu_r', size=1.5)
    pass
