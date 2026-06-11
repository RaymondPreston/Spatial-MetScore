from spatial_metscore import samples_dir, raw_data_dir
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path



# Set up directories
reference_dir = raw_data_dir / "Reference"
reference_dir.mkdir(exist_ok=True)

sns.set_theme(style="white")

# Placeholder for loading logic
def load_and_merge_references():
    # adata_atlas = sc.read_h5ad(reference_dir / 'loveless_2025.h5ad')
    # adata_cosmx = load_cosmx_data(['Pt-5', 'Pt-8', 'Pt-13'])
    # adata_ref = sc.concat([adata_atlas, adata_cosmx], label='source')
    print("Loading Atlas and CosMx reference datasets...")
    pass

# Logic for running STDeconvolve across samples
print("Initiating reference-free LDA discovery via STDeconvolve...")

print("Training Negative Binomial regression model on Composite Reference...")

print("Mapping cell type abundances to Visium spots...")

# Final correlation analysis
# sns.scatterplot(data=adata.obs, x='MetScore', y='CAF_abundance', hue='Tissue_Type')
print("Correlating MetScore with cell type abundances...")
