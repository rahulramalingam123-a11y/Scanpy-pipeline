"""Synthetic count data with known cell types, used for tests and offline demos."""
from __future__ import annotations

import numpy as np
import anndata as ad
import pandas as pd


def make_synthetic_counts(n_cells_per_type: int = 150, n_types: int = 4, n_genes: int = 1200,
                          n_markers: int = 25, n_mito: int = 13, n_bad: int = 20, seed: int = 0):
    """Negative-binomial-like counts with n_types cell types, marker genes, and low-quality cells.

    The first n_mito genes are named MT-*. `n_bad` extra cells have high mitochondrial fraction and
    low depth. Ground-truth labels are in obs['true_type'] ('bad' for the low-quality cells).
    """
    rng = np.random.default_rng(seed)
    base = rng.gamma(0.6, 0.6, size=n_genes)
    base[:n_mito] = 0.3
    blocks, labels = [], []
    for t in range(n_types):
        mu = base.copy()
        marker_idx = n_mito + t * n_markers + np.arange(n_markers)
        mu[marker_idx] *= 12
        depth = rng.lognormal(0, 0.25, size=(n_cells_per_type, 1))
        lam = mu[None, :] * depth
        blocks.append(rng.poisson(rng.gamma(4, lam / 4)))
        labels += [f"type{t}"] * n_cells_per_type
    for _ in range(n_bad):
        mu = base * 0.15
        mu[:n_mito] = base[:n_mito] * 25
        blocks.append(rng.poisson(mu[None, :]))
        labels.append("bad")
    X = np.vstack(blocks).astype(np.float32)
    names = [f"MT-G{i}" if i < n_mito else f"G{i}" for i in range(n_genes)]
    obs = pd.DataFrame({"true_type": labels}, index=[f"cell{i}" for i in range(X.shape[0])])
    adata = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=names))
    return adata
