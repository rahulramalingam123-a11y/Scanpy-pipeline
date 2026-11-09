"""Run scflow on the public 10x PBMC 3k dataset (needs internet; downloads ~6 MB on first run).

    python examples/pbmc3k_demo.py

Outputs go to results_pbmc3k/. Cluster identities are not assigned automatically; use the marker
table (markers.csv) and known PBMC markers (e.g. CD3E for T cells, MS4A1 for B cells, NKG7 for NK
cells, LYZ/CST3 for monocytes) to annotate the clusters yourself.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from scflow import PipelineConfig, run_pipeline
from scflow.pipeline import marker_table

out = Path("results_pbmc3k")
out.mkdir(exist_ok=True)

adata = sc.datasets.pbmc3k()
res = run_pipeline(adata, PipelineConfig(resolution=0.5))
print(res.uns["qc_summary"])
print("clusters:", res.obs["leiden"].nunique())

marker_table(res, 10).to_csv(out / "markers.csv", index=False)
sc.pl.umap(res, color="leiden", show=False)
plt.savefig(out / "umap_leiden.png", dpi=150, bbox_inches="tight")
sc.pl.umap(res, color=["CD3E", "MS4A1", "NKG7", "LYZ"], show=False, ncols=2)
plt.savefig(out / "umap_markers.png", dpi=150, bbox_inches="tight")
