# pgdiffusion (vendored)

Pseudotime Graph Diffusion — random-walk feature diffusion on pseudotime graphs.
Vendored into this repo because it is used **only** inside moma (by
`src/run_pgd_iter.py` and `src/upload_pgd.py`) to build the `X_pseudotime`
embedding behind the PGD UMAP layout.

## Source

- Upstream: https://github.com/brandonlukas/pgdiffusion
- Vendored state: commit `2a0762c` **plus the then-uncommitted working-tree
  version of `diffuse.py`** (the state moma was actually run against).
- Vendored on 2026-07-07.

## Why the working copy, not the release

moma calls `pgd.diffuse(..., return_numpy=True)`. That keyword exists only in the
working-copy `diffuse.py`; the committed release would raise `TypeError`. The
**core diffusion math is identical** between the two — the working copy only adds
input coercion (numpy/list inputs), the `return_numpy` / `show_progress` /
`anchor_mask` / `partition` kwargs, and validation helpers — so results are
numerically unchanged; the working copy is simply the API moma calls.

## Contents

Package only (`__init__.py`, `build_graph.py`, `diffuse.py`). External runtime
deps (torch, numpy, anndata, scanpy, rich) come from the environment.

`import pgdiffusion` resolves here because the scripts run from the moma root
(on `sys.path`); no install step is required.

Note: a separate, later package `pgd2` supersedes this for new work; this vendored
copy is pinned to reproduce the published macrophages analysis.
