"""End-to-end workflow: QC -> filter -> normalise -> HVGs -> PCA -> neighbours -> UMAP -> Leiden -> markers."""
from __future__ import annotations

from dataclasses import dataclass, asdict

import scanpy as sc
from anndata import AnnData

from .qc import add_qc_metrics, filter_cells_genes


@dataclass
class PipelineConfig:
    n_mads: float = 5.0
    max_pct_mt: float | None = 20.0
    min_cells: int = 3
    n_top_genes: int = 2000
    n_pcs: int = 30
    n_neighbors: int = 15
    resolution: float = 0.5
    seed: int = 0
    mito_prefix: str = "MT-"
    n_marker_genes: int = 10


def run_pipeline(adata: AnnData, config: PipelineConfig | None = None) -> AnnData:
    """Run the standard workflow and return a processed copy.

    Raw counts are kept in layers['counts']. Parameters are stored in uns['scflow_config'] so a
    result can be reproduced from the object alone.
    """
    cfg = config or PipelineConfig()
    adata = adata.copy()
    adata.var_names_make_unique()
    add_qc_metrics(adata, cfg.mito_prefix)
    adata = filter_cells_genes(adata, cfg.n_mads, cfg.max_pct_mt, cfg.min_cells)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(cfg.n_top_genes, adata.n_vars - 1))
    adata.raw = adata
    hv = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hv, max_value=10)
    n_pcs = min(cfg.n_pcs, hv.n_vars - 1, hv.n_obs - 1)
    sc.tl.pca(hv, n_comps=n_pcs, random_state=cfg.seed)
    adata.obsm["X_pca"] = hv.obsm["X_pca"]
    sc.pp.neighbors(adata, n_neighbors=cfg.n_neighbors, n_pcs=n_pcs, random_state=cfg.seed)
    sc.tl.umap(adata, random_state=cfg.seed)
    sc.tl.leiden(adata, resolution=cfg.resolution, random_state=cfg.seed,
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)
    adata.uns["scflow_config"] = {k: ("none" if v is None else v) for k, v in asdict(cfg).items()}
    return adata


def marker_table(adata: AnnData, n: int = 10):
    """Top-n marker genes per cluster as a tidy DataFrame (names, scores, log fold changes, adjusted p)."""
    return sc.get.rank_genes_groups_df(adata, group=None).groupby("group", observed=True).head(n)
