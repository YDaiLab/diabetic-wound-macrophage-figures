"""Revision R2-2 - alternative trajectory method (PAGA) cross-check.

The reviewer suggested "comparison with alternative trajectory inference
methods". PAGA (partition-based graph abstraction; Wolf et al. 2019) is an
independent topology method: it scores the connectivity of every cluster pair
from the single-cell kNN graph, with no MST and no pseudotime assumption.

We run PAGA on the JOINT clusters (so results are directly comparable across
groups) for the joint data and, independently, within each condition's cells.
If the joint trajectory reflects real structure, the joint MST edges should be
the high-connectivity cluster pairs under PAGA, in the joint data and within
each condition.

Output:
  results/revision/paga_connectivity.csv
    group (joint/wt/db), source, target, connectivity, in_joint_mst
Console: mean PAGA connectivity on joint-MST edges vs other pairs, and how many
of the joint MST edges fall in each group's top-|MST| PAGA edges.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import kneighbors_graph


def set_knn_graph(a, n_neighbors=15):
    """Attach a symmetric kNN graph (from X_pca) as scanpy neighbors, bypassing
    scanpy's transformer auto-selection (version-flaky here). PAGA operates on
    this graph's connectivities."""
    X = a.obsm["X_pca"]
    conn = kneighbors_graph(X, n_neighbors, mode="connectivity", include_self=False)
    conn = conn.maximum(conn.T)  # symmetric
    dist = kneighbors_graph(X, n_neighbors, mode="distance", include_self=False)
    a.obsp["connectivities"] = conn
    a.obsp["distances"] = dist
    a.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {"n_neighbors": n_neighbors, "method": "umap", "use_rep": "X_pca"},
    }

ADATA_PATH = "box/results/celloracle/cells.h5ad"
EDGE_CSV = "box/results/revision/condpca_edge_detection.csv"  # joint MST edges live here
OUT = "results/revision/paga_connectivity.csv"
SEED = 0
Path(OUT).parent.mkdir(parents=True, exist_ok=True)

# joint MST edges (canonical "a-b", a<b) from the edge-detection table
jedges = pd.read_csv(EDGE_CSV)["edge"].unique().tolist()
jedge_set = {tuple(sorted(map(int, e.split("-")))) for e in jedges}

adata0 = sc.read_h5ad(ADATA_PATH)
adata0 = adata0[~adata0.obs["clusterid"].isna()].copy()
adata0.obs["clusterid"] = adata0.obs["clusterid"].astype(int).astype(str).astype("category")

groups = {
    "joint": adata0,
    "wt": adata0[adata0.obs["condition"] == "wt"].copy(),
    "db": adata0[adata0.obs["condition"] == "db"].copy(),
}

rows = []
for grp, a in groups.items():
    a = a.copy()
    a.obs["clusterid"] = a.obs["clusterid"].cat.remove_unused_categories()
    set_knn_graph(a, n_neighbors=15)
    sc.tl.paga(a, groups="clusterid")

    conn = a.uns["paga"]["connectivities"].toarray()
    cats = list(a.obs["clusterid"].cat.categories)
    for i in range(len(cats)):
        for j in range(i + 1, len(cats)):
            a_, b_ = sorted((int(cats[i]), int(cats[j])))
            rows.append(dict(
                group=grp, source=a_, target=b_,
                connectivity=float(conn[i, j]),
                in_joint_mst=(a_, b_) in jedge_set,
            ))

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(f"wrote {OUT}  ({len(df)} cluster-pairs across {len(groups)} groups)\n")

# summary: are joint MST edges the high-connectivity pairs?
for grp in groups:
    g = df[df.group == grp]
    on = g[g.in_joint_mst].connectivity
    off = g[~g.in_joint_mst].connectivity
    n_mst = int(g.in_joint_mst.sum())
    topN = set(map(tuple, g.nlargest(n_mst, "connectivity")[["source", "target"]].values))
    mst = set(map(tuple, g[g.in_joint_mst][["source", "target"]].values))
    print(f"[{grp:5s}] PAGA connectivity  MST edges mean={on.mean():.3f}  "
          f"other mean={off.mean():.3f}  |  MST edges in top-{n_mst} PAGA: "
          f"{len(topN & mst)}/{n_mst}")
