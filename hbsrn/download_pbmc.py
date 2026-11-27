from pathlib import Path
import scanpy as sc

Path("data").mkdir(exist_ok=True)

adata = sc.datasets.pbmc3k()
adata.write_h5ad("data/pbmc3k_raw.h5ad")

print(adata)
# to read files later
# adata = sc.read_h5ad("data/pbmc3k_raw.h5ad")