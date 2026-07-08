"""Revision R2-2 - per-condition PGD layout.

Runs Pseudotime Graph Diffusion independently within each condition, on that
condition's OWN variant-C trajectory and its OWN PC embedding (nothing borrowed
from the joint model), then UMAPs the diffused features to a condition-specific
PGD layout. Mirrors the production joint run (upload_pgd.py):
neighbors_per_side=50, alpha=0.6, n_steps=1, add_self_loops=True.

Inputs (from export_condpca_for_pgd.R):
  results/revision/pgd_inputs/{wt,db}_orders.csv   branch orders
  results/revision/pgd_inputs/{wt,db}_pca.csv      condition-specific PCA
Output:
  results/revision/pgd_{wt,db}_layout.csv          cell, umap_1, umap_2
"""

import os
import sys
from pathlib import Path

sys.path.insert(
    0, str(Path(__file__).resolve().parent.parent.parent)
)  # moma root -> local pgdiffusion

import numpy as np
import pandas as pd
import scanpy as sc
import pgdiffusion as pgd

ADATA_PATH = "box/results/celloracle/cells.h5ad"
IN_DIR = "results/revision/pgd_inputs"
OUT_DIR = "results/revision"
SEED = 0

adata_all = sc.read_h5ad(ADATA_PATH)

for cond in ["wt", "db"]:
    pca = pd.read_csv(f"{IN_DIR}/{cond}_pca.csv", index_col="cell")
    orders = pd.read_csv(f"{IN_DIR}/{cond}_orders.csv")

    # Subset adata to this condition's cells, in the PCA row order.
    cells = [c for c in pca.index if c in set(adata_all.obs_names)]
    missing = len(pca.index) - len(cells)
    if missing:
        print(f"[{cond}] warning: {missing} trajectory cells not in adata; skipped")
    a = adata_all[cells].copy()
    X = pca.loc[a.obs_names].to_numpy(dtype="float32")

    # branch -> ordered cell ids
    trajectories = {
        b: g.sort_values("pos")["cell"].tolist() for b, g in orders.groupby("branch")
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
    ).assign(condition=cond)
    out.index.name = "cell"
    out_path = f"{OUT_DIR}/pgd_{cond}_layout.csv"
    out.to_csv(out_path)
    print(f"[{cond}] wrote {out_path}  ({a.n_obs} cells, X={X.shape[1]} PCs)")
