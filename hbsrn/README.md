# Single-cell RNA-seq analysis with Scanpy

This script loads an external `.h5ad` file and performs quality control, filtering, normalisation, PCA, Leiden clustering and UMAP visualisation.

The input should contain raw gene-expression counts in `adata.X` or a named layer such as `adata.layers["counts"]`.

## Installation

Use Python 3.12. Open a terminal in the project folder.

### Windows PowerShell

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

### Linux or macOS

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

## Run the analysis

Save the script as `run_external.py` and place your dataset inside a `data` folder.

If raw counts are stored in `adata.X`:

```bash
python run_external.py data/my_dataset.h5ad
```

If raw counts are stored in a layer named `counts`:

```bash
python run_external.py data/my_dataset.h5ad --counts-layer counts
```

Use your environment’s Python executable when running these commands. For example, on Windows:

```powershell
.\.venv\Scripts\python.exe run_external.py data/my_dataset.h5ad
```

To change the filtering threshold and output folder:

```bash
python run_external.py data/my_dataset.h5ad --min-genes 300 --output outputs/sample_01
```

## Analysis steps

1. Load the `.h5ad` file and make gene names unique.
2. Calculate total counts and detected genes for each cell or nucleus.
3. Retain profiles with at least 200 detected genes by default.
4. Retain genes detected in at least three profiles.
5. Preserve the filtered raw counts in `layers["counts"]`.
6. Normalise each profile to 10,000 counts and apply `log1p`.
7. Select up to 2,000 highly variable genes for PCA.
8. Build a neighbour graph, identify Leiden clusters and calculate UMAP coordinates.
9. Save the results.

## Output files

| File | Contents |
|---|---|
| `qc.csv` | QC measurements and filtering decisions for all input profiles |
| `cluster_summary.csv` | Number of retained profiles in each Leiden cluster |
| `umap.png` | UMAP coloured by Leiden cluster |
| `processed.h5ad` | Processed expression, retained raw counts, annotations, clusters and embeddings |

## Interpreting the results

Review the QC distribution before choosing a detected-gene threshold. The default of 200 is a starting value and should be adjusted for the dataset.

Confirm that the selected input matrix contains raw counts before running the script, because it applies normalisation and a log transformation.

Cluster numbers identify groups with similar expression profiles. Assigning biological cell types requires marker-based annotation. Doublet assessment, batch correction and differential-expression analysis are additional steps beyond this script.

## Reference

[Scanpy preprocessing and clustering tutorial](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html)