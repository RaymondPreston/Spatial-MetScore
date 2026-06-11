import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from scipy.stats import rankdata

def  RankAndScore(adata, upSet, downSet, layer="X"):
    """
    Calculate a rank-based gene signature score (singscore) natively in Python.

    This function replicates the logic of the R package 'singscore' (Foroutan et al., 2018).
    It ranks genes within each spot/cell and calculates a normalized score based on 
    provided 'Up' and 'Down' gene signatures. The resulting scores are centered around 0.

    Parameters
    ----------
    adata : AnnData
        Annotated data object containing spatial transcriptomics or scRNA-seq data.
    upSet : list
        List of gene symbols for the 'Up' signature (e.g., Met-high genes).
    downSet : list
        List of gene symbols for the 'Down' signature (e.g., Met-low genes).
    layer : str, optional
        The layer in adata to use for expression values. Defaults to "X".

    Returns
    -------
    AnnData
        The input adata object with 'MetScore', 'MetScore_up', and 'MetScore_down' 
        added to adata.obs.
    
    Notes
    -----
    The score calculation follows these steps:
    1. Rank all genes in a cell by expression level.
    2. Calculate the average rank of genes in the signature.
    3. Normalize the average rank against the theoretical minimum and maximum 
       possible average ranks (S_min and S_max).
    4. Center the scores so that a neutral score is 0.
    """

    # 1. Access the expression data. We convert to a dense array for efficient row-wise ranking.
    # Assumes data in adata.X is already log-normalized.
    expr = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    print(f"Processing matrix of shape: {expr.shape}")

    n_spots, n_total = expr.shape
    gene_names = list(adata.var_names)

    # 2. Filter gene sets to only include those present in the dataset's features.
    met_high_present = [g for g in upSet if g in gene_names]
    met_low_present  = [g for g in downSet  if g in gene_names]

    print(f"Met-high genes: {len(met_high_present)}/{len(upSet)} found in data")
    print(f"Met-low genes:  {len(met_low_present)}/{len(downSet)} found in data")

    # Raise error if a signature is completely missing (prevents division by zero or invalid stats)
    if len(met_high_present) == 0:
        raise ValueError("No met-high genes found in adata.var_names. Check gene symbol format.")
    if len(met_low_present) == 0:
        raise ValueError("No met-low genes found in adata.var_names. Check gene symbol format.")

    # 3. Map gene symbols to their integer indices in the expression matrix.
    gene_idx  = {g: i for i, g in enumerate(gene_names)}
    high_idx  = np.array([gene_idx[g] for g in met_high_present])
    low_idx   = np.array([gene_idx[g] for g in met_low_present])

    n_up   = len(high_idx)
    n_down = len(low_idx)

    # 4. Precompute theoretical boundaries for normalization.
    # S_min: The lowest possible average rank (if signature genes were the lowest expressed).
    # S_max: The highest possible average rank (if signature genes were the highest expressed).
    S_min_up   = (n_up + 1) / 2
    S_max_up   = (2 * n_total - n_up + 1) / 2
    S_min_down = (n_down + 1) / 2
    S_max_down = (2 * n_total - n_down + 1) / 2

    # 5. Rank genes within each spot.
    # We use increasing ranks (highest expression gets highest rank).
    # The down-set uses the exact same ranks, but its score is flipped after normalization.
    ranks = np.apply_along_axis(rankdata, axis=1, arr=expr, method='min')

    # 6. Calculate the observed average rank (S) for the signature genes in each spot.
    S_up   = ranks[:, high_idx].mean(axis=1)
    S_down = ranks[:, low_idx ].mean(axis=1)

    # 7. Normalize the scores to a 0-1 range.
    # Score = (Observed - Min) / (Max - Min)
    up_score   = (S_up   - S_min_up)   / (S_max_up   - S_min_up)
    down_score = (S_down - S_min_down) / (S_max_down - S_min_down)

    # 8. Center and adjust the scores.
    # up_score is shifted by -0.5 so range is [-0.5, 0.5].
    # down_score is flipped (1 - score) and shifted by -0.5 so that high expression of 
    # 'down' genes results in a more negative total score.
    up_score = up_score - 0.5
    down_score = 1 - down_score - 0.5
    
    # Total score range: [-1.0, 1.0]
    total_score = up_score + down_score

    # 9. Store results back into the AnnData object.
    adata.obs['MetScore']      = total_score
    adata.obs['MetScore_up']   = up_score
    adata.obs['MetScore_down'] = down_score

    # Print summary statistics for the user.
    print(f"\nMetScore summary (n={n_spots} spots):")
    print(pd.Series(total_score).describe().round(4).to_string())

    return adata
