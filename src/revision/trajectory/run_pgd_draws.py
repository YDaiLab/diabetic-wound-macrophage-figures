"""Revision R2-2 - per-DRAW PGD layout (10 raw draws x 2 conditions).

Same procedure as run_pgd_condition.py, looped over every raw draw: PGD on the
draw's own variant-C trajectory + its own full PCA, then UMAP.

Inputs (from export_draws_for_pgd.R):
  results/revision/pgd_inputs_draws/{cond}_draw{NN}_{orders,pca}.csv
Output:
  results/revision/pgd_draws/{cond}_draw{NN}_layout.csv   cell, umap_1, umap_2
"""

import sys
from pathlib import Path

sys.path.insert(
    0, str(Path(__file__).resolve().parents[3])
)  # moma root (src/revision/trajectory/ -> moma) -> local pgdiffusion

import pandas as pd
import scanpy as sc
import pgdiffusion as pgd

ADATA_PATH = "box/results/celloracle/cells.h5ad"
IN_DIR = "results/revision/pgd_inputs_draws"
OUT_DIR = "results/revision/pgd_draws"
SEED = 0
Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

adata_all = sc.read_h5ad(ADATA_PATH)
obs_set = set(adata_all.obs_names)

for cond in ["wt", "db"]:
    for d in range(1, 11):
        tag = f"{cond}_draw{d:02d}"
        pca = pd.read_csv(f"{IN_DIR}/{tag}_pca.csv", index_col="cell")
        orders = pd.read_csv(f"{IN_DIR}/{tag}_orders.csv")

        cells = [c for c in pca.index if c in obs_set]
        a = adata_all[cells].copy()
        X = pca.loc[a.obs_names].to_numpy(dtype="float32")
        trajectories = {
            b: g.sort_values("pos")["cell"].tolist()
            for b, g in orders.groupby("branch")
        }

        edge_index = pgd.build_graph(a, trajectories, neighbors_per_side=50)
        X_smooth = pgd.diffuse(
            X, edge_index, alpha=0.4, n_steps=1, add_self_loops=True, return_numpy=True
        )
        a.obsm["X_pseudotime"] = X_smooth
        sc.pp.neighbors(a, use_rep="X_pseudotime", random_state=SEED)
        sc.tl.umap(a, random_state=SEED)

        out = pd.DataFrame(
            a.obsm["X_umap"], columns=["umap_1", "umap_2"], index=a.obs_names
        )
        out.index.name = "cell"
        out.to_csv(f"{OUT_DIR}/{tag}_layout.csv")
        print(f"{tag}: {a.n_obs} cells, X={X.shape[1]} PCs")
