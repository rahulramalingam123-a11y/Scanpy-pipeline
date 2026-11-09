"""Command line: scflow run INPUT --out DIR"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from .pipeline import PipelineConfig, run_pipeline, marker_table
from .simulate import make_synthetic_counts


def _load(path: str):
    p = Path(path)
    if p.is_dir():
        return sc.read_10x_mtx(p, var_names="gene_symbols")
    if p.suffix == ".h5ad":
        return sc.read_h5ad(p)
    if p.suffix == ".h5":
        return sc.read_10x_h5(p)
    raise SystemExit(f"Unsupported input: {path} (use a 10x folder, .h5, or .h5ad)")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="scflow", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    run = sub.add_parser("run", help="Run the workflow on a dataset")
    run.add_argument("input", help="10x mtx folder, 10x .h5, .h5ad, or 'synthetic'")
    run.add_argument("--out", default="scflow_out")
    run.add_argument("--resolution", type=float, default=0.5)
    run.add_argument("--n-mads", type=float, default=5.0)
    run.add_argument("--max-pct-mt", type=float, default=20.0)
    args = ap.parse_args(argv)

    adata = make_synthetic_counts() if args.input == "synthetic" else _load(args.input)
    cfg = PipelineConfig(resolution=args.resolution, n_mads=args.n_mads, max_pct_mt=args.max_pct_mt)
    res = run_pipeline(adata, cfg)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    res.write_h5ad(out / "processed.h5ad")
    marker_table(res, cfg.n_marker_genes).to_csv(out / "markers.csv", index=False)
    sc.pl.umap(res, color="leiden", show=False)
    plt.savefig(out / "umap_leiden.png", dpi=150, bbox_inches="tight")
    plt.close("all")
    print("QC summary:", res.uns["qc_summary"])
    print("Clusters:", res.obs["leiden"].nunique(), "| outputs in", out)


if __name__ == "__main__":
    main()
