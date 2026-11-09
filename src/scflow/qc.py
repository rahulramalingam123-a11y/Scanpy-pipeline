"""Quality control: per-cell metrics and MAD-based outlier filtering."""
from __future__ import annotations

import numpy as np
import scanpy as sc
from anndata import AnnData


def add_qc_metrics(adata: AnnData, mito_prefix: str = "MT-") -> AnnData:
    """Add n_genes_by_counts, total_counts and pct_counts_mt to adata.obs (in place)."""
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(mito_prefix.upper())
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    return adata


def mad_outliers(values, n_mads: float = 5.0, direction: str = "both") -> np.ndarray:
    """Boolean mask of outliers: more than n_mads median absolute deviations from the median.

    direction: 'both', 'high' or 'low'. MAD is scaled by 1.4826 to be comparable to an SD.
    """
    if direction not in ("both", "high", "low"):
        raise ValueError("direction must be 'both', 'high' or 'low'")
    x = np.asarray(values, dtype=float)
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    if mad == 0:
        return np.zeros(x.shape, dtype=bool)
    z = (x - med) / mad
    if direction == "high":
        return z > n_mads
    if direction == "low":
        return z < -n_mads
    return np.abs(z) > n_mads


def filter_cells_genes(adata: AnnData, n_mads: float = 5.0, max_pct_mt: float | None = 20.0,
                       min_cells: int = 3) -> AnnData:
    """Flag and remove low-quality cells, then drop rarely detected genes.

    A cell is removed if log1p(total_counts) or log1p(n_genes) is a low/high MAD outlier, or
    if pct_counts_mt exceeds max_pct_mt (when given). Returns a filtered copy; the number of
    cells removed for each reason is stored in adata.uns['qc_summary'].
    """
    if "pct_counts_mt" not in adata.obs:
        add_qc_metrics(adata)
    out_counts = mad_outliers(np.log1p(adata.obs["total_counts"]), n_mads)
    out_genes = mad_outliers(np.log1p(adata.obs["n_genes_by_counts"]), n_mads)
    out_mt = (adata.obs["pct_counts_mt"].to_numpy() > max_pct_mt) if max_pct_mt is not None \
        else np.zeros(adata.n_obs, dtype=bool)
    remove = out_counts | out_genes | out_mt
    summary = {"n_cells_before": int(adata.n_obs), "outlier_total_counts": int(out_counts.sum()),
               "outlier_n_genes": int(out_genes.sum()), "high_mito": int(out_mt.sum()),
               "removed_any": int(remove.sum())}
    filtered = adata[~remove].copy()
    sc.pp.filter_genes(filtered, min_cells=min_cells)
    summary["n_cells_after"] = int(filtered.n_obs)
    summary["n_genes_after"] = int(filtered.n_vars)
    filtered.uns["qc_summary"] = summary
    return filtered
