# # Tissue-Level MetScore Analysis
# 
# This notebook analyzes the `MetScore` across different tissue types (Primary PDAC, Liver metastasis, Lung metastasis, Peritoneal metastasis) to identify which tissues harbor the most intense metastatic spots.
# 
# We approach this in two ways to account for tumor heterogeneity:
# 1. **Top 95th Percentile Analysis**: Isolating the most aggressive spots in each sample and comparing these peaks across tissue types.
# 2. **Global Distribution (Violin Plots)**: Aggregating every single spot across all samples to visualize the full distribution of scores by tissue type.

from spatial_metscore import samples_dir, metadata_file, data_dir
from spatial_metscore import analyze_sample
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt


output_path = data_dir / "processed"



# 1. Identify Available Samples and Metadata
md_df = pd.read_csv(metadata_file)
available_samples = [d.name for d in samples_dir.iterdir() if d.is_dir() and (d / 'filtered_feature_bc_matrix.h5').exists()]

md_df['Sample_ID_upper'] = md_df['Sample_ID'].str.upper()
md_df = md_df[md_df['Sample_ID_upper'].isin([s.upper() for s in available_samples])]

print(f"Found {len(available_samples)} available samples to process.")

# 2. Extract Scores
# We will build two DataFrames:
# - sample_stats_df: One row per sample (contains 95th percentile score)
# - all_spots_df: One row per spot (contains raw score, for violin plots)

sample_records = []
spot_records = []

for row in md_df.itertuples():
    sample_id_upper = row.Sample_ID_upper
    sample_id = next(s for s in available_samples if s.upper() == sample_id_upper)
    tissue_type = row.Tissue_Type
    patient = row.patient
    sex = row.Sex

    print(f"Processing {sample_id} ({tissue_type})...")
    try:
        adata, _ = analyze_sample(sample_id)
        scores = adata.obs['MetScore'].values

        # Calculate 95th percentile for this sample
        p95_score = np.percentile(scores, 95)

        sample_records.append({
            'Sample_ID': sample_id,
            'Patient': patient,
            'Tissue_Type': tissue_type,
            'Sex': sex,
            'MetScore_95th': p95_score,
            'MetScore_Mean': np.mean(scores),
            'Total_Spots': len(scores)
        })

        # Store all spots (subsampling could be done here if RAM is an issue, but ~150k floats is fine)
        spot_df = pd.DataFrame({
            'Sample_ID': sample_id,
            'Tissue_Type': tissue_type,
            'Sex': sex,
            'MetScore': scores
        })
        spot_records.append(spot_df)

        del adata
        gc.collect()
    except Exception as e:
        print(f"Failed to process {sample_id}: {e}")

sample_stats_df = pd.DataFrame(sample_records)
all_spots_df = pd.concat(spot_records, ignore_index=True)

print(f"\nExtraction complete. Total spots aggregated: {len(all_spots_df)}")


# ### 1. Top Percentile Analysis (The "Hotspot" Average)
# Comparing the 95th percentile of MetScores across tissue types. This shows us the peak aggressiveness achieved in each sample.

# In[ ]:


plt.figure(figsize=(10, 8))
ax = sns.boxplot(data=sample_stats_df, x='Tissue_Type', y='MetScore_95th', palette='Set2', showfliers=False)
sns.stripplot(data=sample_stats_df, x='Tissue_Type', y='MetScore_95th', color='black', alpha=0.6, jitter=True, ax=ax)

plt.title('95th Percentile MetScore by Tissue Type', fontsize=16)
plt.xlabel('Tissue Type', fontsize=12)
plt.ylabel('MetScore (95th Percentile)', fontsize=12)
plt.xticks(rotation=45, ha='right')
sns.despine()

# Wilcoxon rank-sum test and flat line annotation
primary_label = 'Primary PDAC'
unique_types = list(sample_stats_df['Tissue_Type'].unique())
primary_idx = unique_types.index(primary_label)
met_types = [t for t in unique_types if t != primary_label]

primary_data = sample_stats_df[sample_stats_df['Tissue_Type'] == primary_label]['MetScore_95th']
y_max_base = sample_stats_df['MetScore_95th'].max()
line_step = y_max_base * 0.1  # vertical spacing between lines

for i, met in enumerate(met_types):
    met_data = sample_stats_df[sample_stats_df['Tissue_Type'] == met]['MetScore_95th']
    if len(met_data) > 0 and len(primary_data) > 0:
        stat, p_val = mannwhitneyu(primary_data, met_data, alternative='two-sided')

        # Determine asterisks
        if p_val < 0.0001: stars = '****'
        elif p_val < 0.001: stars = '***'
        elif p_val < 0.01: stars = '**'
        elif p_val < 0.05: stars = '*'
        else: stars = 'ns'

        # Calculate positions
        met_idx = unique_types.index(met)
        line_y = y_max_base + (i + 1) * line_step

        # Draw the comparison line (flat, no vertical ticks)
        plt.plot([primary_idx, met_idx], [line_y, line_y], lw=1.5, color='black')

        # Add stars
        plt.text((primary_idx + met_idx) / 2, line_y + line_step*0.05, stars, 
                 ha='center', va='bottom', color='black', fontsize=12)

        p_str = "< 1e-300" if p_val == 0.0 else f"{p_val:.4e}"
        print(f"  Primary PDAC vs {met}: p-value = {p_str} ({stars})")

plt.tight_layout()
plt.savefig(output_path / 'MetScore_95th_Percentile_BoxPlot.png', dpi=300)
plt.show()


# ### 2. Global Distribution (Violin Plots)
# Visualizing the raw distribution of every single spot across all samples, grouped by Tissue Type.

# In[2]:


plt.figure(figsize=(12, 10))
sns.violinplot(data=all_spots_df, x='Tissue_Type', y='MetScore', palette='Set2', inner='quartile', cut=0)

plt.title('Global Distribution of All Spot MetScores by Tissue Type', fontsize=16)
plt.xlabel('Tissue Type', fontsize=12)
plt.ylabel('MetScore', fontsize=12)
plt.xticks(rotation=45, ha='right')
sns.despine()

# Wilcoxon rank-sum test and flat line annotation
primary_label = 'Primary PDAC'
unique_types = list(all_spots_df['Tissue_Type'].unique())
primary_idx = unique_types.index(primary_label)
met_types = [t for t in unique_types if t != primary_label]

primary_data = all_spots_df[all_spots_df['Tissue_Type'] == primary_label]['MetScore']
y_max_base = all_spots_df['MetScore'].max()
line_step = y_max_base * 0.1  # vertical spacing between lines

for i, met in enumerate(met_types):
    met_data = all_spots_df[all_spots_df['Tissue_Type'] == met]['MetScore']
    if len(met_data) > 0 and len(primary_data) > 0:
        stat, p_val = mannwhitneyu(primary_data, met_data, alternative='two-sided')

        # Determine asterisks
        if p_val < 0.0001: stars = '****'
        elif p_val < 0.001: stars = '***'
        elif p_val < 0.01: stars = '**'
        elif p_val < 0.05: stars = '*'
        else: stars = 'ns'

        # Calculate positions
        met_idx = unique_types.index(met)
        line_y = y_max_base + (i + 1) * line_step

        # Draw the comparison line (flat, no vertical ticks)
        plt.plot([primary_idx, met_idx], [line_y, line_y], lw=1.5, color='black')

        # Add stars
        plt.text((primary_idx + met_idx) / 2, line_y + line_step*0.05, stars, 
                 ha='center', va='bottom', color='black', fontsize=12)

        p_str = "< 1e-300" if p_val == 0.0 else f"{p_val:.4e}"
        print(f"  Primary PDAC vs {met}: p-value = {p_str} ({stars})")

plt.tight_layout()
plt.savefig(output_path / 'MetScore_All_Spots_ViolinPlot.png', dpi=300)
plt.show()


# In[ ]:




