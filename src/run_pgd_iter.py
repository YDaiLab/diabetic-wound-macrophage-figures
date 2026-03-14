import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
import scanpy as sc
import pgdiffusion as pgd
from utils import load_orders_from_lamian
import scipy.sparse as sp

adata = sc.read_h5ad("box/results/celloracle/cells.h5ad")
trajectories = load_orders_from_lamian("box/results/lamian/infer_tree.rds")

df_0 = pd.DataFrame(
    adata.obsm["X_umap"], columns=["umap_1", "umap_2"], index=adata.obs_names
).assign(alpha=0)

# Get embeddings (PCA, UMAP, etc.)
X = adata.obsm["X_pca"]

# Build pseudotime graph
edge_index = pgd.build_graph(adata, trajectories, neighbors_per_side=50)

# Apply diffusion
X_smooth = pgd.diffuse(
    X,
    edge_index,
    alpha=0.2,
    n_steps=1,
    add_self_loops=True,
    return_numpy=True,
)

# Store results
adata.obsm["X_pseudotime"] = X_smooth

# Visualize
sc.pp.neighbors(adata, use_rep="X_pseudotime")
sc.tl.umap(adata)

df_2 = pd.DataFrame(
    adata.obsm["X_umap"], columns=["umap_1", "umap_2"], index=adata.obs_names
).assign(alpha=0.2)

# Apply diffusion
X_smooth = pgd.diffuse(
    X,
    edge_index,
    alpha=0.4,
    n_steps=1,
    add_self_loops=True,
    return_numpy=True,
)

# Store results
adata.obsm["X_pseudotime"] = X_smooth

# Visualize
sc.pp.neighbors(adata, use_rep="X_pseudotime")
sc.tl.umap(adata)

df_4 = pd.DataFrame(
    adata.obsm["X_umap"], columns=["umap_1", "umap_2"], index=adata.obs_names
).assign(alpha=0.4)

# Load other adata
adata = sc.read_h5ad("box/results/pgd_celloracle/cells.h5ad")
df_6 = pd.DataFrame(
    adata.obsm["X_umap"], columns=["umap_1", "umap_2"], index=adata.obs_names
).assign(alpha=0.6)

# Combine all dataframes
df_all = pd.concat([df_0, df_2, df_4, df_6]).merge(
    adata.obs, left_index=True, right_index=True, how="left", suffixes=("", "_obs")
)
df_all.to_csv("results/pgd_iter_meta.csv")
