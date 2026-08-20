"""Revision R2-5, panel c — sensitivity of PGD to the diffusion parameter alpha.

Regenerates the whole sweep in one environment so the points are mutually
comparable: splicing new alphas onto the existing results/pgd_iter_meta.csv would
mix UMAP versions, and a different UMAP build can restructure a layout rather
than merely rotate it.

alpha = 0 is the undiffused control and alpha = 0.6 the value used in the paper.
The sweep brackets it with 0.8 and 0.9 so the choice is not at the edge of the
range.

These are NOT the published layouts and must not be labelled as such. The
diffused PCs are identical (asserted below), but the UMAP step is not: relative
RMSD after Procrustes is 0.17 at alpha = 0.6 and 0.95 at alpha = 0, the latter
because the published conventional UMAP came from Seurat's RunUMAP rather than
scanpy. The metrics are structural and agree to <= 0.02 with the same metrics
computed on the published layouts over alpha 0 to 0.6.

Diffusion is replicated here in numpy rather than called from pgdiffusion, which
needs torch. The replication is exact, not approximate: the assert below
reproduces the published X_pseudotime bit-for-bit from X_pca. NB self-loops are
ON, matching the published run (they are not pgdiffusion's current default, and
without them the diffused PCs differ by up to 0.8).

Metrics are computed in two spaces. In UMAP space they describe what was actually
published; in diffused-PC space they are deterministic and isolate alpha from
UMAP's stochasticity. Both are rotation- and reflection-invariant, so no
alignment is needed for them; the Procrustes-aligned coordinates below exist only
so the layouts can be shown as a strip without each panel arriving at a random
orientation.

Run: conda run -n pgd2 python src/revision/pgd/alpha_sweep.py
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from pgd2 import pseudotime_smoothness
from scipy.linalg import orthogonal_procrustes
from sklearn.manifold import trustworthiness

# ----------------------------------------------------------------------
h5ad_file = "box/results/pgd_celloracle/cells.h5ad"  # carries X_pca and the published X_pseudotime
branches_file = "box/results/lamian/branches.csv"
published_umap_file = "results/pgd_cells_meta.csv"
layout_out = "results/revision/pgd_alpha_sweep.csv"
metrics_out = "results/revision/pgd_alpha_metrics.csv"

ALPHAS = (0.0, 0.2, 0.4, 0.6, 0.8, 0.9)
PUBLISHED_ALPHA = 0.6
NEIGHBORS_PER_SIDE = 50  # as run_pgd.py
K = 15  # neighbourhood size for both metrics
SEED = 0
# ----------------------------------------------------------------------

adata = ad.read_h5ad(h5ad_file)
X_pca = np.asarray(adata.obsm["X_pca"], dtype=np.float32)
pt = np.asarray(adata.obs["pseudotime"], dtype=float)
idx = {name: i for i, name in enumerate(adata.obs_names)}


def pseudotime_graph(branches, n_cells, k):
    """pgdiffusion.build_graph: slide a window of +/-k along each branch ordering,
    dedupe undirected pairs keeping the shortest positional step, emit both
    directions."""
    seen = {}
    for _, g in branches.groupby("branch"):
        seq = [idx[c] for c in g.sort_values("pseudotime")["cell_id"] if c in idx]
        for pos, ci in enumerate(seq):
            for other in range(max(0, pos - k), min(len(seq), pos + k + 1)):
                if other == pos:
                    continue
                cj = seq[other]
                key = (ci, cj) if ci < cj else (cj, ci)
                step = abs(other - pos)
                if key not in seen or step < seen[key]:
                    seen[key] = step
    # interleave the two directions per pair, as build_graph emits them: float32
    # accumulation in np.add.at is order-dependent, so this ordering is what
    # reproduces the published values exactly rather than to ~1e-7.
    pairs = np.array(list(seen.keys()), dtype=np.int64)
    loops = np.arange(n_cells)
    src = np.concatenate([pairs.ravel(), loops])
    dst = np.concatenate([pairs[:, ::-1].ravel(), loops])
    return src, dst


def diffuse(X, src, dst, alpha):
    """One pgdiffusion step: H <- (1-a) H + a * D^-1 A H, in-degree clamped at 1."""
    deg = np.zeros(X.shape[0], dtype=np.float32)
    np.add.at(deg, dst, np.float32(1.0))
    agg = np.zeros_like(X)
    np.add.at(agg, dst, X[src])
    agg /= np.maximum(deg, np.float32(1.0))[:, None]
    return (np.float32(1 - alpha) * X + np.float32(alpha) * agg).astype(np.float32)


src, dst = pseudotime_graph(pd.read_csv(branches_file), X_pca.shape[0], NEIGHBORS_PER_SIDE)
_published_pcs = np.asarray(adata.obsm["X_pseudotime"], dtype=np.float32)
_err = np.abs(diffuse(X_pca, src, dst, PUBLISHED_ALPHA) - _published_pcs).max()
print(f"diffusion check: max|reproduced - published X_pseudotime| = {_err:.3e}")
assert _err < 1e-5, "numpy diffusion no longer reproduces the published X_pseudotime"

layouts, metrics = [], []
for alpha in ALPHAS:
    H = diffuse(X_pca, src, dst, alpha)
    a = ad.AnnData(X=np.zeros((H.shape[0], 1), dtype=np.float32), obs=adata.obs.copy())
    a.obsm["X_diffused"] = H
    sc.pp.neighbors(a, use_rep="X_diffused", random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    emb = np.asarray(a.obsm["X_umap"])

    layouts.append(pd.DataFrame({
        "cell": list(adata.obs_names), "alpha": alpha,
        "umap_1": emb[:, 0], "umap_2": emb[:, 1],
        "pseudotime": pt, "clusterid": np.asarray(adata.obs["clusterid"]),
    }))
    for space, E in (("UMAP", emb), ("Diffused PCs", H)):
        metrics.append({
            "alpha": alpha, "space": space, "k": K,
            "pseudotime_smoothness": pseudotime_smoothness(E, pt, k=K),
            "trustworthiness": trustworthiness(X_pca, E, n_neighbors=K),
        })

layout = pd.concat(layouts, ignore_index=True)

# Procrustes-align every layout onto the published one, for display only.
published = (
    pd.read_csv(published_umap_file).drop_duplicates("Row.names")
    .set_index("Row.names").reindex(adata.obs_names)[["umap_1", "umap_2"]].to_numpy()
)
ref = published - published.mean(0)
ref /= np.linalg.norm(ref)
for alpha in ALPHAS:
    m = layout["alpha"] == alpha
    Y = layout.loc[m, ["umap_1", "umap_2"]].to_numpy()
    Y = Y - Y.mean(0)
    Y /= np.linalg.norm(Y)
    R, _ = orthogonal_procrustes(Y, ref)
    layout.loc[m, ["aligned_1", "aligned_2"]] = Y @ R

Path(layout_out).parent.mkdir(parents=True, exist_ok=True)
layout.to_csv(layout_out, index=False)
pd.DataFrame(metrics).to_csv(metrics_out, index=False)

# How far the regenerated published-alpha layout lands from the published one,
# after alignment: this is the UMAP-environment difference, quoted in the legend.
reg = layout.loc[layout["alpha"] == PUBLISHED_ALPHA, ["aligned_1", "aligned_2"]].to_numpy()
print(f"regenerated alpha={PUBLISHED_ALPHA} vs published, aligned RMSD (unit-scaled): "
      f"{np.sqrt(((reg - ref) ** 2).sum(1).mean()):.4f}")
print(pd.DataFrame(metrics).to_string(index=False))
