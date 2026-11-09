import numpy as np
import pytest
from sklearn.metrics import adjusted_rand_score

from scflow import (PipelineConfig, add_qc_metrics, filter_cells_genes, make_synthetic_counts,
                    mad_outliers, run_pipeline)
from scflow.pipeline import marker_table
from scflow.cli import main


def test_mad_outliers_flags_extremes_only():
    x = np.r_[np.random.default_rng(0).normal(10, 1, 200), [100.0, -80.0]]
    m = mad_outliers(x, 5)
    assert m[-2:].all() and m[:-2].sum() == 0


def test_mad_direction_and_constant_input():
    x = np.r_[np.ones(50), [50.0]]
    assert not mad_outliers(np.ones(10)).any()
    assert mad_outliers(np.r_[np.random.default_rng(1).normal(0, 1, 100), [50.0]], 5, "high")[-1]
    assert not mad_outliers(np.r_[np.random.default_rng(1).normal(0, 1, 100), [50.0]], 5, "low")[-1]
    with pytest.raises(ValueError):
        mad_outliers(x, 5, "sideways")


def test_qc_metrics_mito_fraction():
    a = make_synthetic_counts(n_bad=10, seed=1)
    add_qc_metrics(a)
    bad = a.obs["true_type"] == "bad"
    assert a.obs.loc[bad, "pct_counts_mt"].mean() > 3 * a.obs.loc[~bad, "pct_counts_mt"].mean()


def test_filter_removes_bad_cells_and_keeps_good():
    a = make_synthetic_counts(n_bad=25, seed=2)
    f = filter_cells_genes(a, n_mads=5, max_pct_mt=20)
    kept_bad = (f.obs["true_type"] == "bad").sum()
    kept_good = (f.obs["true_type"] != "bad").sum()
    assert kept_bad <= 3
    assert kept_good >= 0.95 * (a.obs["true_type"] != "bad").sum()
    s = f.uns["qc_summary"]
    assert s["n_cells_after"] == f.n_obs and s["removed_any"] == a.n_obs - f.n_obs


def test_pipeline_recovers_known_cell_types():
    a = make_synthetic_counts(seed=3)
    r = run_pipeline(a, PipelineConfig(n_top_genes=500, resolution=0.5))
    truth = r.obs["true_type"]
    assert (truth == "bad").sum() <= 3
    ari = adjusted_rand_score(truth[truth != "bad"], r.obs.loc[truth != "bad", "leiden"])
    assert ari > 0.9
    assert "counts" in r.layers and "X_umap" in r.obsm
    assert r.uns["scflow_config"]["resolution"] == 0.5


def test_markers_are_the_planted_genes():
    a = make_synthetic_counts(n_markers=25, seed=4)
    r = run_pipeline(a, PipelineConfig(n_top_genes=500))
    mk = marker_table(r, 10)
    # every top marker should be a planted marker gene (index >= 13 and < 13+4*25)
    idx = mk["names"].str.replace("G", "", regex=False).astype(int)
    assert ((idx >= 13) & (idx < 13 + 4 * 25)).mean() > 0.9


def test_pipeline_is_deterministic():
    a = make_synthetic_counts(seed=5)
    r1 = run_pipeline(a, PipelineConfig(n_top_genes=500))
    r2 = run_pipeline(a, PipelineConfig(n_top_genes=500))
    assert (r1.obs["leiden"].to_numpy() == r2.obs["leiden"].to_numpy()).all()


def test_cli_synthetic(tmp_path):
    main(["run", "synthetic", "--out", str(tmp_path / "o")])
    assert (tmp_path / "o" / "processed.h5ad").exists()
    assert (tmp_path / "o" / "markers.csv").exists()
    assert (tmp_path / "o" / "umap_leiden.png").exists()
