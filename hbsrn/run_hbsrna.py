"""Run QC and clustering on an external h5ad count matrix."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path, help="Input h5ad file")
    parser.add_argument("--counts-layer", help="Layer containing raw counts")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--output", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    output = args.output
    output.mkdir(parents=True, exist_ok=True)

    # Load your external dataset.
    raw = sc.read_h5ad(args.input)
    raw.var_names_make_unique()

    # Use this option if counts are stored in a layer.
    if args.counts_layer:
        raw.X = raw.layers[args.counts_layer].copy()

    print(f"Loaded {raw.n_obs} cells/nuclei × {raw.n_vars} genes")

    sc.pp.calculate_qc_metrics(raw, percent_top=None, inplace=True)
    keep = raw.obs["n_genes_by_counts"] >= args.min_genes
    raw.obs.assign(passed_qc=keep).to_csv(output / "qc.csv")

    data = raw[keep].copy()
    sc.pp.filter_genes(data, min_cells=3)

    if min(data.shape) < 3:
        raise ValueError("Too few cells or genes remain after filtering.")

    data.layers["counts"] = data.X.copy()
    sc.pp.normalize_total(data, target_sum=10_000)
    sc.pp.log1p(data)

    # Select informative genes for PCA.
    sc.pp.highly_variable_genes(
        data,
        n_top_genes=min(2_000, data.n_vars),
        flavor="seurat",
    )

    n_variable = int(data.var["highly_variable"].sum())
    n_pcs = min(30, data.n_obs - 1, n_variable - 1)

    if n_pcs < 2:
        raise ValueError("Too few variable genes remain for PCA.")

    sc.pp.pca(data, n_comps=n_pcs, random_state=42)
    sc.pp.neighbors(
        data,
        n_neighbors=min(12, data.n_obs - 1),
        n_pcs=n_pcs,
        random_state=42,
    )
    sc.tl.leiden(
        data,
        flavor="igraph",
        directed=False,
        n_iterations=2,
        resolution=0.4,
        random_state=42,
    )
    sc.tl.umap(data, random_state=42)

    sc.pl.umap(data, color="leiden", show=False)
    plt.savefig(output / "umap.png", dpi=160, bbox_inches="tight")
    plt.close("all")

    data.obs["leiden"].value_counts().sort_index().rename_axis(
        "cluster"
    ).to_csv(output / "cluster_summary.csv", header=["n_cells"])

    data.write_h5ad(output / "processed.h5ad")
    print(f"Retained {data.n_obs} profiles. Results saved to {output.resolve()}")


if __name__ == "__main__":
    main()