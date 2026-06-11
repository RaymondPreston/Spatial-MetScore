# # MetScore Analysis & Spatial Visualization
#
# This script processes all available Visium spatial transcriptomics samples to:
#   1. Calculate MetScore for each spot (with per-sample .h5ad caching)
#   2. Generate side-by-side H&E / MetScore spatial plots grouped by patient
#   3. Compare MetScore distributions across tissue types (boxplot + violin)
#
# Tissue types: Primary PDAC, Liver metastasis, Lung metastasis, Peritoneal metastasis

import gc
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import mannwhitneyu

from spatial_metscore import samples_dir, metadata_file, data_dir, analyze_sample

# ── Paths ─────────────────────────────────────────────────────────────────────
processed_dir = data_dir / "processed"
cache_dir     = processed_dir / "h5ad_cache"
figures_dir   = processed_dir / "figures"

processed_dir.mkdir(parents=True, exist_ok=True)
cache_dir.mkdir(parents=True, exist_ok=True)
figures_dir.mkdir(parents=True, exist_ok=True)

# ── Constants ─────────────────────────────────────────────────────────────────
MAX_ROWS_PER_SLIDE = 3   # max samples per patient figure
METSCORE_VMIN      = -0.1
METSCORE_VMAX      =  0.3


# ══════════════════════════════════════════════════════════════════════════════
# Identify Available Samples
# ══════════════════════════════════════════════════════════════════════════════

md_df = pd.read_csv(metadata_file)
available_samples = [
    d.name for d in samples_dir.iterdir()
    if d.is_dir() and (d / "filtered_feature_bc_matrix.h5").exists()
]

md_df["Sample_ID_upper"] = md_df["Sample_ID"].str.upper()
md_df = md_df[md_df["Sample_ID_upper"].isin([s.upper() for s in available_samples])]

patients = md_df["patient"].unique()
print(f"Found {len(available_samples)} samples across {len(patients)} patients.")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Load or Process Each Sample (with .h5ad cache)
# ══════════════════════════════════════════════════════════════════════════════

def load_or_process(sample_id, load_images=False):
    """
    Return (adata, metadata) for a sample.
    Loads from .h5ad cache if available; otherwise runs analyze_sample() and caches.
    Images are never cached — pass load_images=True only when needed for spatial plots.
    """
    cache_path = cache_dir / f"{sample_id}.h5ad"

    if cache_path.exists() and not load_images:
        print(f"  [cache] Loading {sample_id} from {cache_path.name}")
        adata = sc.read_h5ad(cache_path)
        # Retrieve metadata from the master table
        match = md_df[md_df["Sample_ID_upper"] == sample_id.upper()]
        metadata = match.iloc[0].to_dict() if not match.empty else {}
        return adata, metadata

    print(f"  [process] Running analyze_sample for {sample_id}...")
    adata, metadata = analyze_sample(sample_id, load_images=load_images)

    if not load_images:
        adata.write_h5ad(cache_path)
        print(f"  [cache] Saved to {cache_path.name}")

    return adata, metadata


# ══════════════════════════════════════════════════════════════════════════════
# 3. Extract Per-Sample and Per-Spot Statistics
# ══════════════════════════════════════════════════════════════════════════════

sample_records = []
spot_records   = []

for row in md_df.itertuples():
    sample_id    = next(s for s in available_samples if s.upper() == row.Sample_ID_upper)
    tissue_type  = row.Tissue_Type
    patient      = row.patient
    sex          = row.Sex

    print(f"Processing {sample_id} ({tissue_type})...")
    try:
        adata, _ = load_or_process(sample_id, load_images=False)
        scores   = adata.obs["MetScore"].values

        sample_records.append({
            "Sample_ID":    sample_id,
            "Patient":      patient,
            "Tissue_Type":  tissue_type,
            "Sex":          sex,
            "MetScore_95th": np.percentile(scores, 95),
            "MetScore_Mean": np.mean(scores),
            "Total_Spots":   len(scores),
        })

        spot_records.append(pd.DataFrame({
            "Sample_ID":   sample_id,
            "Tissue_Type": tissue_type,
            "Sex":         sex,
            "MetScore":    scores,
        }))

        del adata
        gc.collect()

    except Exception as e:
        print(f"  Failed to process {sample_id}: {e}")

sample_stats_df = pd.DataFrame(sample_records)
all_spots_df    = pd.concat(spot_records, ignore_index=True)

print(f"\nExtraction complete. Total spots: {len(all_spots_df)}")

# Save aggregated tables
sample_stats_df.to_csv(processed_dir / "sample_stats.csv", index=False)
all_spots_df.to_parquet(processed_dir / "all_spots.parquet")


# ══════════════════════════════════════════════════════════════════════════════
# 4. Statistical Annotation Helper
# ══════════════════════════════════════════════════════════════════════════════

def annotate_pairwise(ax, data_df, group_col, value_col, reference_group):
    """
    Draw Wilcoxon rank-sum significance bars comparing each group to reference_group.

    Parameters
    ----------
    ax             : matplotlib Axes
    data_df        : DataFrame containing group_col and value_col
    group_col      : column name for the categorical grouping variable
    value_col      : column name for the numeric values
    reference_group: the group to compare all others against (e.g. 'Primary PDAC')
    """
    unique_types = list(data_df[group_col].unique())
    ref_idx      = unique_types.index(reference_group)
    ref_data     = data_df[data_df[group_col] == reference_group][value_col]

    y_max_base = data_df[value_col].max()
    line_step  = y_max_base * 0.1

    comparison_count = 0
    for other_group in unique_types:
        if other_group == reference_group:
            continue

        other_data = data_df[data_df[group_col] == other_group][value_col]
        if len(other_data) == 0 or len(ref_data) == 0:
            continue

        stat, p_val = mannwhitneyu(ref_data, other_data, alternative="two-sided")

        if   p_val < 0.0001: stars = "****"
        elif p_val < 0.001:  stars = "***"
        elif p_val < 0.01:   stars = "**"
        elif p_val < 0.05:   stars = "*"
        else:                stars = "ns"

        other_idx = unique_types.index(other_group)
        line_y    = y_max_base + (comparison_count + 1) * line_step

        ax.plot([ref_idx, other_idx], [line_y, line_y], lw=1.5, color="black")
        ax.text(
            (ref_idx + other_idx) / 2, line_y + line_step * 0.05,
            stars, ha="center", va="bottom", color="black", fontsize=12,
        )

        p_str = "< 1e-300" if p_val == 0.0 else f"{p_val:.4e}"
        print(f"  {reference_group} vs {other_group}: p = {p_str} ({stars})")
        comparison_count += 1


# ══════════════════════════════════════════════════════════════════════════════
# 5. Plot A — 95th Percentile Boxplot (one point per sample)
# ══════════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(10, 8))
sns.boxplot(
    data=sample_stats_df, x="Tissue_Type", y="MetScore_95th",
    palette="Set2", showfliers=False, ax=ax,
)
sns.stripplot(
    data=sample_stats_df, x="Tissue_Type", y="MetScore_95th",
    color="black", alpha=0.6, jitter=True, ax=ax,
)
annotate_pairwise(ax, sample_stats_df, "Tissue_Type", "MetScore_95th", "Primary PDAC")

ax.set_title("95th Percentile MetScore by Tissue Type", fontsize=16)
ax.set_xlabel("Tissue Type", fontsize=12)
ax.set_ylabel("MetScore (95th Percentile)", fontsize=12)
ax.tick_params(axis="x", rotation=45)
sns.despine()
plt.tight_layout()
plt.savefig(figures_dir / "MetScore_95th_Percentile_BoxPlot.png", dpi=300)
plt.show()


# ══════════════════════════════════════════════════════════════════════════════
# 6. Plot B — Global Distribution Violin (all spots)
# ══════════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(12, 10))
sns.violinplot(
    data=all_spots_df, x="Tissue_Type", y="MetScore",
    palette="Set2", inner="quartile", cut=0, ax=ax,
)
annotate_pairwise(ax, all_spots_df, "Tissue_Type", "MetScore", "Primary PDAC")

ax.set_title("Global Distribution of All Spot MetScores by Tissue Type", fontsize=16)
ax.set_xlabel("Tissue Type", fontsize=12)
ax.set_ylabel("MetScore", fontsize=12)
ax.tick_params(axis="x", rotation=45)
sns.despine()
plt.tight_layout()
plt.savefig(figures_dir / "MetScore_All_Spots_ViolinPlot.png", dpi=300)
plt.show()


# ══════════════════════════════════════════════════════════════════════════════
# 7. Spatial Visualization — H&E + MetScore per Patient
# ══════════════════════════════════════════════════════════════════════════════

def generate_patient_slide(patient_id, chunk_index, samples_chunk):
    """
    Generate and save a figure with H&E (left) and MetScore overlay (right)
    for up to MAX_ROWS_PER_SLIDE samples of a given patient.
    """
    n_rows = len(samples_chunk)
    fig, axs = plt.subplots(n_rows, 2, figsize=(16, 6 * n_rows))
    axs = np.atleast_2d(axs)   # always 2D, safe for n_rows == 1

    for i, row_data in enumerate(samples_chunk.itertuples()):
        sample_id = next(s for s in available_samples if s.upper() == row_data.Sample_ID_upper)
        print(f"  Rendering {sample_id}...")

        try:
            adata, metadata = load_or_process(sample_id, load_images=True)

            tissue_type = metadata.get("Tissue_Type", "Unknown Tissue")
            treatment   = metadata.get("treatment", "") or ""
            title_prefix = f"{sample_id} ({patient_id}) - {tissue_type}\n{treatment}".strip()

            # Left: H&E
            sq.pl.spatial_scatter(
                adata, color=None, ax=axs[i, 0],
                title=f"{title_prefix}\nH&E Stain",
                na_color=(0, 0, 0, 0),
            )

            # Right: MetScore overlay
            sq.pl.spatial_scatter(
                adata, color="MetScore", cmap="RdBu_r",
                size=1.5, vmin=METSCORE_VMIN, vmax=METSCORE_VMAX,
                ax=axs[i, 1],
                title=f"{title_prefix}\nMetScore",
                colorbar=True,
            )

        except Exception as e:
            print(f"  Failed to render {sample_id}: {e}")
            axs[i, 0].set_title(f"Failed: {sample_id}\n{e}")
            axs[i, 1].set_title("Error")

    plt.tight_layout()

    suffix   = f"_part{chunk_index}" if chunk_index > 1 else ""
    out_path = figures_dir / f"{patient_id}{suffix}.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out_path.name}\n")
    return out_path


generated_plots = []

for patient in patients:
    print(f"====== {patient} ======")
    patient_samples = md_df[md_df["patient"] == patient]
    n_chunks = (len(patient_samples) + MAX_ROWS_PER_SLIDE - 1) // MAX_ROWS_PER_SLIDE

    for chunk_idx in range(n_chunks):
        chunk = patient_samples.iloc[chunk_idx * MAX_ROWS_PER_SLIDE : (chunk_idx + 1) * MAX_ROWS_PER_SLIDE]
        plot_path = generate_patient_slide(patient, chunk_idx + 1, chunk)
        generated_plots.append(plot_path)

print(f"\nDone. {len(generated_plots)} spatial figures saved to {figures_dir}")
