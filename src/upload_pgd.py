import scanpy as sc
import pgdiffusion as pgd
from utils import load_orders_from_lamian
from pathlib import Path

adata = sc.read_h5ad("box/results/celloracle/cells.h5ad")
trajectories = load_orders_from_lamian("box/results/lamian/infer_tree.rds")

# Get embeddings (PCA, UMAP, etc.)
X = adata.obsm["X_pca"]

# Build pseudotime graph
edge_index = pgd.build_graph(adata, trajectories, neighbors_per_side=50)

# Apply diffusion
X_smooth = pgd.diffuse(
    X,
    edge_index,
    alpha=0.6,
    n_steps=1,
    add_self_loops=True,
    return_numpy=True,
)

# Store results
adata.obsm["X_pseudotime"] = X_smooth

# Visualize
sc.pp.neighbors(adata, use_rep="X_pseudotime")
sc.tl.umap(adata)
sc.pl.umap(adata, color=["pseudotime", "clusterid"])

# output_file = Path(
#     "/Users/brandonlukas/Library/CloudStorage/Box-Box/data/timkoh/macrophages/from_mac/pgd_cells.h5ad"
# )
# output_file.parent.mkdir(parents=True, exist_ok=True)
# adata.write_h5ad(output_file)
