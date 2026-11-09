"""Offline demo on synthetic data with known cell types; no download needed.

    python examples/synthetic_demo.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

from scflow import make_synthetic_counts, run_pipeline

adata = make_synthetic_counts(seed=0)
res = run_pipeline(adata)
ok = res.obs["true_type"] != "bad"
print(res.uns["qc_summary"])
print("clusters:", res.obs["leiden"].nunique())
print("adjusted Rand index vs truth:", round(adjusted_rand_score(res.obs.loc[ok, "true_type"], res.obs.loc[ok, "leiden"]), 3))
fig, ax = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"wspace": 0.45})
sc.pl.umap(res, color="leiden", ax=ax[0], show=False, title="Leiden clusters")
sc.pl.umap(res, color="true_type", ax=ax[1], show=False, title="True cell type")
plt.savefig("docs/synthetic_demo.png", dpi=150, bbox_inches="tight")
