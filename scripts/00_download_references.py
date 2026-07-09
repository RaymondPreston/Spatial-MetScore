"""
Reference Dataset Downloader
=============================
Downloads processed, annotated scRNA-seq files for building the Tangram
composite reference. Run from the project root with the Spatial conda env active.

Usage:
    python scripts/00_download_references.py

Datasets:
    1. Loveless 2025  — PDAC atlas          (Zenodo 14199536,  ~33 GB)
    2. Peng 2019      — PDAC primary        (Zenodo 3969339,   ~2.8 GB)
    3. Raghavan 2021  — Metastatic PDAC     (Broad SCP1644,    MANUAL)
    4. HLCA           — Human Lung Atlas    (CellxGene CDN,    ~21 GB)
    5. Lambrechts 2018— NSCLC TME           (HCA / ArrayExpress, ~loom files)

Requirements:
    pip install zenodo-get requests tqdm
"""

import os
import sys
import subprocess
import textwrap
import requests
from pathlib import Path
from tqdm import tqdm

# ── Paths ─────────────────────────────────────────────────────────────────────
# Resolve project root relative to this script's location
SCRIPT_DIR   = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
REF_DIR      = PROJECT_ROOT / "data" / "raw" / "Reference"
REF_DIR.mkdir(parents=True, exist_ok=True)


# ══════════════════════════════════════════════════════════════════════════════
# Utility helpers
# ══════════════════════════════════════════════════════════════════════════════

def wget_download(url: str, dest: Path, label: str = "") -> bool:
    """Download a file with wget, showing progress. Skip if already exists."""
    if dest.exists():
        print(f"  [skip] {dest.name} already exists.")
        return True
    print(f"  [download] {label or dest.name}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        ["wget", "-c", "-O", str(dest), url],
        check=False,
    )
    if result.returncode != 0:
        print(f"  [error] wget failed for {url}")
        return False
    print(f"  [done] Saved to {dest.name}")
    return True


def zenodo_download(record_id: str, target_dir: Path) -> bool:
    """Download all files from a Zenodo record using zenodo_get."""
    target_dir.mkdir(parents=True, exist_ok=True)
    print(f"  [zenodo_get] Downloading record {record_id} → {target_dir}")
    try:
        subprocess.run(
            ["zenodo_get", str(record_id)],
            cwd=str(target_dir),
            check=True,
        )
        return True
    except FileNotFoundError:
        print("  [error] zenodo_get not found. Install with: pip install zenodo-get")
        return False
    except subprocess.CalledProcessError as e:
        print(f"  [error] zenodo_get failed: {e}")
        return False


# ══════════════════════════════════════════════════════════════════════════════
# 1. Loveless 2025 — PDAC Atlas (Zenodo 14199536)
# ══════════════════════════════════════════════════════════════════════════════

def download_loveless():
    """
    Loveless IM et al. (2025) Human pancreatic cancer single cell atlas.
    Clin Cancer Res. Zenodo record 14199536.
    Single h5ad file (~33 GB).
    """
    print("\n══ 1. Loveless 2025 PDAC Atlas ══")
    out_dir = REF_DIR / "loveless_2025"

    # Check if already downloaded (zenodo_get names the file from the record)
    existing = list(out_dir.glob("*.h5ad")) + list(out_dir.glob("*.rds")) + list(out_dir.glob("*.RDS"))
    if existing:
        print(f"  [skip] Found existing file(s): {[f.name for f in existing]}")
        return True

    success = zenodo_download("14199536", out_dir)

    # Standardize filename if needed
    if success:
        rds_files = list(out_dir.glob("*.rds")) + list(out_dir.glob("*.RDS"))
        h5ad_files = list(out_dir.glob("*.h5ad"))
        if h5ad_files:
            target = out_dir / "loveless_2025_atlas.h5ad"
            if h5ad_files[0] != target:
                h5ad_files[0].rename(target)
                print(f"  [rename] → {target.name}")
        elif rds_files:
            print("  [note] Downloaded as .rds — run reference_preparation.py to convert to h5ad")
    return success


# ══════════════════════════════════════════════════════════════════════════════
# 2. Peng 2019 — PDAC Primary (Zenodo 3969339)
# ══════════════════════════════════════════════════════════════════════════════

def download_peng():
    """
    Peng J et al. (2019) Single-cell RNA-seq highlights intra-tumoral
    heterogeneity in PDAC. Cell Research.
    Reprocessed h5ad files from Roche/Besca (Zenodo 3969339, ~2.8 GB total).
    Contains 57,530 cells: ductal, acinar, endocrine, stellate, fibroblast,
    endothelial, macrophage, T cell, B cell.
    """
    print("\n══ 2. Peng 2019 PDAC ══")
    out_dir = REF_DIR / "peng_2019"

    existing = list(out_dir.glob("*.h5ad"))
    if len(existing) >= 2:
        print(f"  [skip] Found {len(existing)} h5ad files already.")
        return True

    # Fetch file list from Zenodo API
    print("  [api] Fetching file list from Zenodo record 3969339...")
    try:
        resp = requests.get(
            "https://zenodo.org/api/records/3969339",
            timeout=30,
        )
        resp.raise_for_status()
        record = resp.json()
        files = record.get("files", [])
        if not files:
            # v2 API format
            files = record.get("metadata", {}).get("files", [])
    except Exception as e:
        print(f"  [error] Could not fetch Zenodo API: {e}")
        print("  [fallback] Using zenodo_get CLI...")
        return zenodo_download("3969339", out_dir)

    out_dir.mkdir(parents=True, exist_ok=True)
    success = True
    for f in files:
        fname = f.get("key") or f.get("filename", "unknown")
        url   = f.get("links", {}).get("self") or f.get("download_url", "")
        if not url:
            continue
        dest = out_dir / fname
        ok = wget_download(url, dest, label=fname)
        success = success and ok

    if not success:
        print("  [fallback] Trying zenodo_get CLI...")
        return zenodo_download("3969339", out_dir)

    return success


# ══════════════════════════════════════════════════════════════════════════════
# 3. Raghavan 2021 — Metastatic PDAC (Broad SCP1644) — MANUAL
# ══════════════════════════════════════════════════════════════════════════════

def download_raghavan_instructions():
    """
    Raghavan S et al. (2021) Microenvironment drives cell state, plasticity,
    and drug response in pancreatic cancer. Cell.
    23 metastatic PDAC biopsies (mostly liver mets). Seq-Well platform.
    Requires a free Broad Single Cell Portal account to download.
    """
    print("\n══ 3. Raghavan 2021 Metastatic PDAC ══")
    out_dir = REF_DIR / "raghavan_2021"
    out_dir.mkdir(parents=True, exist_ok=True)

    instructions_file = out_dir / "README_manual_download.txt"
    instructions = textwrap.dedent("""
    Raghavan 2021 — Manual Download Instructions
    =============================================
    URL: https://singlecell.broadinstitute.org/single_cell/study/SCP1644

    The Broad Single Cell Portal requires a Google account login to download files.
    Automated download is not possible without authentication tokens.

    Steps:
    1. Go to: https://singlecell.broadinstitute.org/single_cell/study/SCP1644
    2. Click "Sign In" (top right) and log in with a Google account.
    3. Click the "Download" tab.
    4. Download the following files into this directory:
       - The main expression matrix (DGE matrix, .txt.gz or similar)
       - The cell metadata file (cell_metadata.txt or similar)
       - The cluster/annotation file if separate

    Alternatively, use the SCP command-line tool:
        pip install scp-cli
        scp-cli auth login
        scp-cli download --study SCP1644 --output-dir .

    After downloading, place all files in:
        data/raw/Reference/raghavan_2021/

    Expected files:
        - expression matrix (genes x cells, sparse or dense)
        - cell metadata with cell type annotations
        - ~23 samples, ~68,000 cells total
    """).strip()

    instructions_file.write_text(instructions)
    print(f"  [manual required] Instructions saved to {instructions_file}")
    print("  ⚠  Raghavan 2021 requires a Broad SCP login — see README in raghavan_2021/")
    return False  # Signals that manual action is needed


# ══════════════════════════════════════════════════════════════════════════════
# 4. HLCA — Human Lung Cell Atlas (CellxGene)
# ══════════════════════════════════════════════════════════════════════════════

def download_hlca():
    """
    Sikkema L et al. (2023) An integrated cell atlas of the human lung in
    health and disease. Nature Medicine.
    2.28M cells, comprehensive lung cell type annotation.
    Direct download from CellxGene CDN (~21 GB).
    """
    print("\n══ 4. HLCA — Human Lung Cell Atlas ══")
    out_dir = REF_DIR / "hlca"
    dest    = out_dir / "hlca_core.h5ad"

    # Confirmed URL from CellxGene API (collection 6f6d381a, dataset dbb5ad81)
    url = "https://datasets.cellxgene.cziscience.com/688185ad-11c2-4172-a53a-f4f1f4076860.h5ad"

    return wget_download(url, dest, label="HLCA core (~21 GB)")


# ══════════════════════════════════════════════════════════════════════════════
# 5. Lambrechts 2018 — NSCLC TME (HCA / ArrayExpress E-MTAB-6149 + E-MTAB-6653)
# ══════════════════════════════════════════════════════════════════════════════

def download_lambrechts():
    """
    Lambrechts D et al. (2018) Phenotype molding of stromal cells in the lung
    tumor microenvironment. Nature Medicine.
    ~93k cells, 52 stromal subtypes, NSCLC TME.
    Loom files available via HCA Data Explorer.
    """
    print("\n══ 5. Lambrechts 2018 NSCLC TME ══")
    out_dir = REF_DIR / "lambrechts_2018"
    out_dir.mkdir(parents=True, exist_ok=True)

    existing_looms = list(out_dir.glob("*.loom*"))
    if existing_looms:
        print(f"  [skip] Found {len(existing_looms)} loom file(s) already.")
        return True

    # Fetch file manifest from HCA DCP API
    print("  [api] Fetching file manifest from HCA Data Explorer...")
    project_id = "453d7ee2-319f-496c-9862-99d397870b63"
    manifest_url = (
        f"https://service.azul.data.humancellatlas.org/index/projects/"
        f"{project_id}?catalog=dcp3"
    )

    try:
        resp = requests.get(manifest_url, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        # Extract loom file URLs from the project manifest
        loom_urls = []
        for file_entry in data.get("files", {}).get("hits", []):
            for f in file_entry.get("files", []):
                if f.get("format", "").lower() in ("loom", "loom.gz"):
                    url = f.get("url", "")
                    name = f.get("name", f"lambrechts_{len(loom_urls)}.loom.gz")
                    if url:
                        loom_urls.append((name, url))

        if loom_urls:
            print(f"  [found] {len(loom_urls)} loom file(s) via HCA API")
            success = True
            for name, url in loom_urls:
                dest = out_dir / name
                ok = wget_download(url, dest, label=name)
                success = success and ok
            return success

    except Exception as e:
        print(f"  [api error] {e}")

    # Fallback: print manual instructions for ArrayExpress
    print("  [fallback] HCA API did not return direct loom URLs.")
    instructions_file = out_dir / "README_manual_download.txt"
    instructions = textwrap.dedent("""
    Lambrechts 2018 — Manual Download Instructions
    ===============================================
    Loom and RDS files are available via the HCA Data Explorer:
    https://explore.data.humancellatlas.org/projects/453d7ee2-319f-496c-9862-99d397870b63

    Steps:
    1. Go to the URL above and click "Download" or "Matrices"
    2. Select the loom.gz files (8 files, ~few GB total)
    3. Place all downloaded files in:
       data/raw/Reference/lambrechts_2018/

    Alternatively, download directly from ArrayExpress:
    - E-MTAB-6149: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6149
    - E-MTAB-6653: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6653
    Look for the processed loom or RDS files in the supplementary files section.

    Note: The loom files contain raw counts + cell type annotations.
    You will need to convert them to h5ad in the harmonization script.
    """).strip()
    instructions_file.write_text(instructions)
    print(f"  [manual required] Instructions saved to {instructions_file}")
    return False


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("Reference Dataset Downloader")
    print(f"Output directory: {REF_DIR}")
    print("=" * 60)

    results = {
        "Loveless 2025 (PDAC atlas)":    download_loveless(),
        "Peng 2019 (PDAC primary)":      download_peng(),
        "Raghavan 2021 (met PDAC)":      download_raghavan_instructions(),
        "HLCA (Human Lung Atlas)":       download_hlca(),
        "Lambrechts 2018 (NSCLC TME)":  download_lambrechts(),
    }

    print("\n" + "=" * 60)
    print("Download Summary")
    print("=" * 60)
    for dataset, ok in results.items():
        status = "✓ done" if ok else "⚠  manual action required"
        print(f"  {status:30s}  {dataset}")

    manual_needed = [k for k, v in results.items() if not v]
    if manual_needed:
        print(f"\n{len(manual_needed)} dataset(s) require manual download.")
        print("Check the README_manual_download.txt files in each dataset folder.")
    else:
        print("\nAll datasets downloaded successfully.")


if __name__ == "__main__":
    main()
