#!/usr/bin/env python
# coding: utf-8

# # MetScore Spatial Visualization Report
# 
# This notebook processes all available spatial transcriptomics samples, calculates the native Python `MetScore`, and generates high-resolution spatial plots.
# 
# Plots are generated side-by-side (H&E on the left, MetScore on the right) and grouped by patient. Each patient gets up to 3 samples (rows) per figure/slide. Plots are saved to the `MetScore Visualization/` directory.

# In[1]:


import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
import os
import gc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from spatial_metscore.config import samples_dir, metadata_file, data_dir
from spatial_metscore import analyze_sample

output_dir = Path('figures')
output_dir.mkdir(exist_ok=True)


# In[10]:


# 1. Group samples by Patient

if metadata_file.exists():
    md_df = pd.read_csv(metadata_file)
else:
    raise FileNotFoundError(f"Cannot find {metadata_file}")

# Find all fully downloaded samples (those that have the h5 matrix)
available_samples = [d.name for d in samples_dir.iterdir() if d.is_dir() and (d / 'filtered_feature_bc_matrix.h5').exists()]

# Filter metadata to only include available samples
# Note: metadata Sample_ID might have different casing than the folder, let's normalize for matching if needed.
# In the metadata it's 'Pt-1A', folder is 'PT-1A'
md_df['Sample_ID_upper'] = md_df['Sample_ID'].str.upper()
md_df = md_df[md_df['Sample_ID_upper'].isin([s.upper() for s in available_samples])]

# Group by patient
patients = md_df['patient'].unique()
print(f"Found {len(available_samples)} available samples across {len(patients)} patients.")


# In[3]:


def generate_patient_slide(patient_id, chunk_index, samples_chunk):
    """
    Generates a single figure (slide) for up to 3 samples of a given patient.
    """
    n_rows = len(samples_chunk)
    # Figure size: width=16 (2 columns), height=6 per row
    fig, axs = plt.subplots(n_rows, 2, figsize=(16, 6 * n_rows))

    # Handle 1D array of axes if n_rows == 1
    if n_rows == 1:
        axs = [axs]

    for i, row_data in enumerate(samples_chunk.itertuples()):
        # Convert from metadata ID (Pt-1A) to folder name (PT-1A)
        sample_id_upper = row_data.Sample_ID_upper
        # Find exact folder name
        sample_id = next(s for s in available_samples if s.upper() == sample_id_upper)

        print(f"Processing {sample_id}...")
        try:
            adata, metadata = analyze_sample(sample_id, load_images=True)

            tissue_type = metadata.get('Tissue_Type', 'Unknown Tissue')
            treatment = metadata.get('treatment', '')
            if pd.isna(treatment):
                treatment = ''

            title_prefix = f"{sample_id} ({patient_id}) - {tissue_type}\n{treatment}"

            # Left: H&E Stain
            sq.pl.spatial_scatter(
                adata, 
                color=None,
                ax=axs[i][0],
                title=f"{title_prefix}\nH&E Stain",
                na_color=(0,0,0,0)  # Transparent spots so just image shows
            )

            # Right: MetScore Overlay
            sq.pl.spatial_scatter(
                adata, 
                color='MetScore', 
                cmap='RdBu_r', 
                size=1.5, 
                vmin=-0.1,
                vmax=0.3,
                ax=axs[i][1],
                title=f"{title_prefix}\nMetScore",
                colorbar=True
            )

            del adata
            gc.collect()
        except Exception as e:
            print(f"Failed to process {sample_id}: {e}")
            axs[i][0].set_title(f"Failed to load {sample_id}\n{e}")
            axs[i][1].set_title("Error")

    plt.tight_layout()

    # Save the plot
    suffix = f"_part{chunk_index}" if chunk_index > 1 else ""
    out_path = output_dir / f"{patient_id}{suffix}.png"
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved slide to {out_path}\n")

    return out_path


# In[ ]:


from IPython.display import Image, display

generated_plots = []
MAX_ROWS = 3

for patient in patients:
    print(f"====== Generating slides for {patient} ======")
    patient_samples = md_df[md_df['patient'] == patient]

    # Chunk into groups of MAX_ROWS (3)
    num_chunks = (len(patient_samples) + MAX_ROWS - 1) // MAX_ROWS

    for chunk_idx in range(num_chunks):
        start_idx = chunk_idx * MAX_ROWS
        end_idx = start_idx + MAX_ROWS
        chunk = patient_samples.iloc[start_idx:end_idx]

        plot_path = generate_patient_slide(patient, chunk_idx + 1, chunk)
        generated_plots.append(plot_path)
        display(Image(filename=plot_path))

