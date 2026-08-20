"""Revision R2-2 - spring (force-directed) layout from BIRDccNEST's own graph.

Panel h shows the BIRDccNEST communities + oMST trajectory. Rather than borrow
the Seurat/PGD UMAP, we lay the cells out with a force-directed layout computed
from BIRDccNEST's OWN pruned cell-cell graph (pruned_adj.npz), so the layout,
communities, flow and backbone all come from the same graph.

The pruned graph is dense (avg degree ~1.5k); a spring layout on it is a
hairball. For the layout only we keep each cell's TOP_K strongest edges
(symmetrized union) - standard sparsification for force-directed single-cell
layouts - then run scanpy's ForceAtlas2 draw_graph (spreads graph structure into
trajectory-like arms). Edges still come entirely from pruned_adj; TOP_K only
controls visual spread.

Inputs:
  box/results/revision/birdccnest/joint_res1p5/{pruned_adj.npz, cell_names.npy}
Output:
  results/revision/birdccnest_spring_layout.csv   cell, umap_1, umap_2
"""

import random
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp

BASE = "box/results/revision/birdccnest/joint_res1p5"
OUT = "results/revision/birdccnest_spring_layout.csv"
EDGES_OUT = "results/revision/birdccnest_spring_edges.csv"
LAYOUT = "fa"  # ForceAtlas2 (fa2_modified); "fr" = Fruchterman-Reingold
TOP_K = 15
SEED = 0
Path(OUT).parent.mkdir(parents=True, exist_ok=True)

random.seed(SEED)  # python-igraph (draw_graph backend) draws from python's RNG
np.random.seed(SEED)


def topk_rows(A, k):
    """Keep the k largest-weight entries per row of a CSR matrix."""
    A = A.tocsr()
    keep_d, keep_i, keep_p = [], [], [0]
    for r in range(A.shape[0]):
        s, e = A.indptr[r], A.indptr[r + 1]
        idx, dat = A.indices[s:e], A.data[s:e]
        if len(dat) > k:
            top = np.argpartition(dat, -k)[-k:]
            idx, dat = idx[top], dat[top]
        keep_i.extend(idx)
        keep_d.extend(dat)
        keep_p.append(keep_p[-1] + len(dat))
    return sp.csr_matrix((keep_d, keep_i, keep_p), shape=A.shape)


A = sp.load_npz(f"{BASE}/pruned_adj.npz")
names = np.load(f"{BASE}/cell_names.npy", allow_pickle=True)

# strongest TOP_K per cell, then symmetric union -> undirected weighted graph
S = topk_rows(A, TOP_K)
S = S.maximum(S.T)
print(f"sparsified {A.nnz} -> {S.nnz} nnz (avg degree {S.nnz / S.shape[0]:.1f})")

adata = ad.AnnData(sp.csr_matrix((A.shape[0], 1), dtype="float32"))
adata.obs_names = [str(n) for n in names]
adata.obsp["connectivities"] = S
adata.obsp["distances"] = S
adata.uns["neighbors"] = {
    "connectivities_key": "connectivities",
    "distances_key": "distances",
    "params": {"n_neighbors": TOP_K, "method": "umap", "use_rep": "birdccnest"},
}

sc.tl.draw_graph(adata, layout=LAYOUT, random_state=SEED)
xy = adata.obsm[f"X_draw_graph_{LAYOUT}"]

import pandas as pd

out = pd.DataFrame({"cell": adata.obs_names, "umap_1": xy[:, 0], "umap_2": xy[:, 1]})
out.to_csv(OUT, index=False)
print(f"wrote {OUT}  ({len(out)} cells)")

# undirected cell-cell edges (upper triangle) with layout coords, for the
# cell-graph panel. Endpoints stored directly so R just draws segments.
St = sp.triu(S, k=1).tocoo()
edges = pd.DataFrame({
    "x": xy[St.row, 0], "y": xy[St.row, 1],
    "xend": xy[St.col, 0], "yend": xy[St.col, 1],
    "w": St.data,
})
edges.to_csv(EDGES_OUT, index=False)
print(f"wrote {EDGES_OUT}  ({len(edges)} edges)")
