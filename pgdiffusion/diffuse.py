import math
from typing import Union
import torch
import numpy as np
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)


def _coerce_tensor(
    value: Union[torch.Tensor, np.ndarray, list, None],
    device: torch.device,
    dtype: torch.dtype,
    name: str,
) -> Union[torch.Tensor, None]:
    """
    Convert input to torch.Tensor on the specified device and dtype.

    Parameters
    ----------
    value : torch.Tensor, np.ndarray, list, or None
        Input value to coerce.
    device : torch.device
        Target device.
    dtype : torch.dtype
        Target dtype.
    name : str
        Parameter name for error messages.

    Returns
    -------
    torch.Tensor or None
        Coerced tensor, or None if input is None.
    """
    if value is None:
        return None
    if isinstance(value, torch.Tensor):
        return value.to(device=device, dtype=dtype)
    if isinstance(value, np.ndarray):
        return torch.from_numpy(value).to(device=device, dtype=dtype)
    if isinstance(value, (list, tuple)):
        return torch.tensor(value, device=device, dtype=dtype)
    raise TypeError(
        f"{name} must be a torch.Tensor, numpy array, or list, got {type(value)}"
    )


def _coerce_tensor_int(
    value: Union[torch.Tensor, np.ndarray, list, None],
    device: torch.device,
    name: str,
) -> Union[torch.Tensor, None]:
    """
    Convert input to torch.Tensor with int64 dtype on the specified device.

    Parameters
    ----------
    value : torch.Tensor, np.ndarray, list, or None
        Input value to coerce.
    device : torch.device
        Target device.
    name : str
        Parameter name for error messages.

    Returns
    -------
    torch.Tensor or None
        Coerced tensor, or None if input is None.
    """
    if value is None:
        return None
    if isinstance(value, torch.Tensor):
        return value.to(device=device, dtype=torch.long)
    if isinstance(value, np.ndarray):
        return torch.from_numpy(value).to(device=device, dtype=torch.long)
    if isinstance(value, (list, tuple)):
        return torch.tensor(value, device=device, dtype=torch.long)
    raise TypeError(
        f"{name} must be a torch.Tensor, numpy array, or list, got {type(value)}"
    )


def _coerce_tensor_bool(
    value: Union[torch.Tensor, np.ndarray, list, None],
    device: torch.device,
    name: str,
) -> Union[torch.Tensor, None]:
    """
    Convert input to torch.Tensor with bool dtype on the specified device.

    Parameters
    ----------
    value : torch.Tensor, np.ndarray, list, or None
        Input value to coerce.
    device : torch.device
        Target device.
    name : str
        Parameter name for error messages.

    Returns
    -------
    torch.Tensor or None
        Coerced tensor, or None if input is None.
    """
    if value is None:
        return None
    if isinstance(value, torch.Tensor):
        return value.to(device=device, dtype=torch.bool)
    if isinstance(value, np.ndarray):
        return torch.from_numpy(value).to(device=device, dtype=torch.bool)
    if isinstance(value, (list, tuple)):
        return torch.tensor(value, device=device, dtype=torch.bool)
    raise TypeError(
        f"{name} must be a torch.Tensor, numpy array, or list, got {type(value)}"
    )


def _validate_diffuse_inputs(
    X: torch.Tensor,
    edge_index: torch.Tensor,
    alpha: float,
    n_steps: int,
    edge_weight: Union[torch.Tensor, None],
    add_self_loops: bool,
    self_loop_weight: float,
    M: Union[torch.Tensor, None],
    U: Union[torch.Tensor, None],
    beta: float,
    anchor_mask: Union[torch.Tensor, None],
    partition: Union[torch.Tensor, None],
) -> None:
    """
    Validate all diffuse() inputs for consistency and correctness.

    Parameters
    ----------
    X : torch.Tensor
        Feature matrix (N, d)
    edge_index : torch.Tensor
        Edge index (2, E)
    alpha : float
        Blend weight
    n_steps : int
        Number of steps
    edge_weight : torch.Tensor or None
        Edge weights
    add_self_loops : bool
        Whether to add self-loops
    self_loop_weight : float
        Self-loop weight
    M : torch.Tensor or None
        Feature mixing matrix
    U : torch.Tensor or None
        Low-rank coupling matrix
    beta : float
        Coupling strength
    anchor_mask : torch.Tensor or None
        Anchor node mask
    partition : torch.Tensor or None
        Node partition labels

    Raises
    ------
    ValueError
        If any input is invalid.
    TypeError
        If input type is invalid.
    """
    N, d = X.shape

    if edge_index.shape[0] != 2:
        raise ValueError(f"edge_index must have shape (2, E), got {edge_index.shape}")
    if edge_index.min() < 0 or edge_index.max() >= N:
        raise ValueError(f"edge_index contains node ids outside [0, {N-1}]")
    if not 0.0 <= alpha <= 1.0:
        raise ValueError(f"alpha must be in [0, 1], got {alpha}")
    if not isinstance(n_steps, int) or n_steps < 1:
        raise ValueError(f"n_steps must be a positive integer, got {n_steps}")
    if add_self_loops:
        if not math.isfinite(self_loop_weight):
            raise ValueError("self_loop_weight must be finite")
        if self_loop_weight < 0:
            raise ValueError("self_loop_weight must be >= 0")
    if M is not None and U is not None:
        raise ValueError("Provide either M (explicit) or U (low-rank), not both.")
    if M is not None:
        if M.dim() != 2 or M.shape[0] != d or M.shape[1] != d:
            raise ValueError(
                f"M must have shape (n_features, n_features) = ({d}, {d}), "
                f"got {tuple(M.shape)}"
            )
    if U is not None:
        if U.dim() != 2 or U.shape[0] != d:
            raise ValueError(f"U must have shape (d, r)=({d}, r), got {tuple(U.shape)}")
        if beta < 0:
            raise ValueError("beta must be >= 0")
    if anchor_mask is not None:
        if anchor_mask.dtype != torch.bool or anchor_mask.shape != (N,):
            raise ValueError(
                f"anchor_mask must be a boolean tensor of shape (N,), got "
                f"dtype={anchor_mask.dtype}, shape={tuple(anchor_mask.shape)}"
            )
    if partition is not None:
        if partition.dim() != 1 or partition.shape[0] != N:
            raise ValueError(
                f"partition must have shape (N,) = ({N},), got {tuple(partition.shape)}"
            )
    if edge_weight is not None:
        if edge_weight.dim() != 1 or edge_weight.shape[0] != edge_index.shape[1]:
            raise ValueError(
                f"edge_weight must be 1D with length equal to number of edges ({edge_index.shape[1]}), "
                f"got shape {tuple(edge_weight.shape)}"
            )
        if not torch.isfinite(edge_weight).all():
            raise ValueError("edge_weight must contain only finite values")
        if (edge_weight < 0).any():
            raise ValueError("edge_weight must be >= 0")


def _prepare_edges(
    edge_index: torch.Tensor,
    edge_weight: Union[torch.Tensor, None],
    add_self_loops: bool,
    self_loop_weight: float,
    partition: Union[torch.Tensor, None],
    N: int,
    device: torch.device,
    dtype: torch.dtype,
) -> tuple:
    """
    Prepare edge lists and weights, including self-loops and partition filtering.

    Self-loops are added before partition filtering.

    Parameters
    ----------
    edge_index : torch.Tensor
        Edge index (2, E)
    edge_weight : torch.Tensor or None
        Edge weights (E,)
    add_self_loops : bool
        Whether to add self-loops
    self_loop_weight : float
        Weight for self-loops
    partition : torch.Tensor or None
        Node partition labels
    N : int
        Number of nodes
    device : torch.device
        Target device
    dtype : torch.dtype
        Target dtype for weights

    Returns
    -------
    tuple
        (src, dst, w) where src and dst are 1D tensors and w is 1D weight tensor
    """
    src, dst = edge_index

    # Prepare edge weights (default to ones)
    if edge_weight is not None:
        w = edge_weight.to(device=device, dtype=dtype)
    else:
        w = torch.ones(src.shape[0], device=device, dtype=dtype)

    # Add self-loops before partition filtering
    if add_self_loops:
        self_edges = torch.arange(N, device=device)
        src = torch.cat([src, self_edges])
        dst = torch.cat([dst, self_edges])
        w = torch.cat(
            [
                w,
                torch.full(
                    (N,),
                    float(self_loop_weight),
                    device=device,
                    dtype=dtype,
                ),
            ]
        )

    # Filter edges by partition
    if partition is not None:
        same = partition[src] == partition[dst]
        src = src[same]
        dst = dst[same]
        w = w[same]

    return src, dst, w


def diffuse(
    X: Union[torch.Tensor, np.ndarray, list],
    edge_index: Union[torch.Tensor, np.ndarray, list],
    *,
    alpha: float = 0.6,
    n_steps: int = 1,
    edge_weight: Union[torch.Tensor, np.ndarray, list, None] = None,
    add_self_loops: bool = False,
    self_loop_weight: float = 1.0,
    M: Union[torch.Tensor, np.ndarray, None] = None,
    U: Union[torch.Tensor, np.ndarray, None] = None,
    beta: float = 0.0,
    anchor_mask: Union[torch.Tensor, np.ndarray, list, None] = None,
    partition: Union[torch.Tensor, np.ndarray, list, None] = None,
    return_numpy: bool = False,
    show_progress: bool = False,
):
    """
    Random-walk feature diffusion on a pseudotime graph.

    We update features by blending a residual term with an incoming-neighbor
    (src -> dst) weighted mean aggregation, followed by optional feature coupling.

    **Mathematical formulation**: In the simplest case (unweighted graph, no feature
    coupling), let H^(t) in R^(N x d) denote the feature matrix at iteration t. One
    diffusion step is given by

        H^(t+1) = (1 - alpha) H^(t) + alpha * P H^(t),

    where the propagation operator P is defined as

        P = D^(-1) A.

    Here A is the adjacency matrix with entries A[i, s] = 1 if there is an edge
    s -> i and 0 otherwise, and D is the diagonal in-degree matrix with
    D[i, i] = sum_s A[i, s].

    Parameters
    ----------
    X : torch.Tensor, numpy.ndarray, or list
        Feature matrix (N, d). Can be a numpy array, list, or torch tensor.
        Internally converted to torch.Tensor.
    edge_index : torch.Tensor, numpy.ndarray, or list
        Sparse edges (2, E), with (src, dst). Can be a numpy array, list, or torch tensor.
        Internally converted to torch.Tensor with dtype=torch.long.
    alpha : float
        Blend weight in [0, 1]. Default: 0.6.
    n_steps : int
        Number of diffusion iterations. Default: 1.
    edge_weight : torch.Tensor, numpy.ndarray, list, or None
        Optional per-edge weights aligned with edge_index columns (E,).
        Can be a numpy array, list, or torch tensor. Internally converted to torch.Tensor.
        Default: None (uniform weights).
    add_self_loops : bool
        If True, include self-loops in the aggregation operator. Default: False.
        Self-loops are added before partition filtering.
    self_loop_weight : float
        Weight assigned to self-loop edges when add_self_loops=True. Default: 1.0.
    M : torch.Tensor, numpy.ndarray, or None
        Optional explicit feature mixing matrix (d, d). Default: None.
    U : torch.Tensor, numpy.ndarray, or None
        Optional low-rank feature coupling (d, r) (e.g., PCA loadings). Default: None.
    beta : float
        Coupling strength for low-rank U; ignored if U is None. Default: 0.0.
    anchor_mask : torch.Tensor, numpy.ndarray, list, or None
        Boolean mask of shape (N,) indicating nodes clamped to original X after
        each iteration. Default: None. Anchored nodes override diffusion at every step.
    partition : torch.Tensor, numpy.ndarray, list, or None
        Integer labels of length N restricting diffusion to within-group edges.
        Edges are filtered after self-loops are added, keeping only edges where
        partition[src] == partition[dst]. Default: None.
    return_numpy : bool
        If True, return H.cpu().numpy(). If False (default), return torch.Tensor.
        Default: False.
    show_progress : bool
        If True, display a progress bar during diffusion. If False, run silently.
        Default: False.

    Returns
    -------
    torch.Tensor or numpy.ndarray
        Smoothed features with the same shape as X.
        Type depends on return_numpy parameter.

    Raises
    ------
    ValueError
        If both M and U are provided, or if shapes are invalid.
    TypeError
        If input types are invalid.

    Examples
    --------
    >>> import torch
    >>> import scanpy as sc
    >>> import pgdiffusion as pgd
    >>> adata = sc.read_h5ad("data.h5ad")
    >>> trajectories = {"branch": ["cell_0", "cell_1", "cell_2"]}
    >>> X = torch.tensor(adata.obsm["X_pca"], dtype=torch.float32)
    >>> edge_index = pgd.build_graph(adata, trajectories)
    >>> edge_index = torch.tensor(edge_index, dtype=torch.long)
    >>> X_smooth = pgd.diffuse(X, edge_index, alpha=0.6, n_steps=1)
    >>> adata.obsm["X_pseudotime"] = X_smooth.cpu().numpy()

    Numpy array input example
    -------------------------
    Input arrays are automatically converted to torch tensors.

    >>> X = np.random.randn(100, 50)
    >>> edge_index = np.array([[0, 1, 2], [1, 2, 3]])
    >>> X_smooth = pgd.diffuse(X, edge_index, alpha=0.6, n_steps=1, return_numpy=True)

    Weighted edges example
    ----------------------
    Use positional steps from `build_graph(..., include_step_attr=True)` to derive
    per-edge weights (e.g., inverse-distance) and pass them via `edge_weight`.

    >>> edge_index, edge_steps = pgd.build_graph(
    ...     adata, trajectories, include_step_attr=True
    ... )
    >>> edge_index = torch.tensor(edge_index, dtype=torch.long)
    >>> edge_weight = torch.tensor(1.0 / (edge_steps + 1.0), dtype=torch.float32)
    >>> X_smooth = pgd.diffuse(
    ...     X,
    ...     edge_index,
    ...     alpha=0.6,
    ...     n_steps=1,
    ...     edge_weight=edge_weight,
    ... )

    Feature coupling (PCA loadings) example
    ---------------------------------------
    Apply low-rank feature coupling using PCA loadings U (e.g., Scanpy stores gene
    loadings in ``adata.varm["PCs"]``). This biases diffusion toward the principal
    subspace while preserving full feature dimensionality.

    >>> X = torch.tensor(adata.X.toarray(), dtype=torch.float32)  # (n_cells, n_genes)
    >>> U = torch.tensor(adata.varm["PCs"], dtype=torch.float32)  # (n_genes, n_pcs)
    >>> X_smooth = pgd.diffuse(X, edge_index, alpha=0.5, n_steps=5, U=U, beta=0.3)

    Anchor mask example
    -------------------
    Keep certain nodes fixed to their original values:

    >>> anchor_mask = torch.zeros(100, dtype=torch.bool)
    >>> anchor_mask[[0, 10, 50]] = True  # Fix nodes 0, 10, 50
    >>> X_smooth = pgd.diffuse(X, edge_index, anchor_mask=anchor_mask, n_steps=5)

    Partition example
    ------------------
    Restrict diffusion to within-group edges (e.g., different cell types):

    >>> partition = torch.tensor([0, 0, 0, 1, 1, 1])  # Two groups
    >>> X_smooth = pgd.diffuse(X, edge_index, partition=partition, n_steps=5)
    """
    # ===== Input coercion =====
    # Coerce X first to determine device and dtype for other inputs
    X_tensor = _coerce_tensor(X, torch.device("cpu"), torch.float32, "X")
    if X_tensor is None:
        raise TypeError("X cannot be None")

    # Coerce other inputs to match X's device and dtype
    device = X_tensor.device
    dtype = X_tensor.dtype

    edge_index_tensor = _coerce_tensor_int(edge_index, device, "edge_index")
    if edge_index_tensor is None:
        raise TypeError("edge_index cannot be None")

    edge_weight_tensor = _coerce_tensor(edge_weight, device, dtype, "edge_weight")
    M_tensor = _coerce_tensor(M, device, dtype, "M")
    U_tensor = _coerce_tensor(U, device, dtype, "U")
    anchor_mask_tensor = _coerce_tensor_bool(anchor_mask, device, "anchor_mask")
    partition_tensor = _coerce_tensor_int(partition, device, "partition")

    # ===== Validation =====
    _validate_diffuse_inputs(
        X_tensor,
        edge_index_tensor,
        alpha,
        n_steps,
        edge_weight_tensor,
        add_self_loops,
        self_loop_weight,
        M_tensor,
        U_tensor,
        beta,
        anchor_mask_tensor,
        partition_tensor,
    )

    N, d = X_tensor.shape

    # ===== Edge preparation =====
    src, dst, w = _prepare_edges(
        edge_index_tensor,
        edge_weight_tensor,
        add_self_loops,
        self_loop_weight,
        partition_tensor,
        N,
        device,
        dtype,
    )

    # ===== Precompute weighted in-degree =====
    deg = torch.zeros(N, device=device, dtype=dtype)
    deg.index_add_(0, dst, w)
    deg = deg.clamp_min(1.0)

    # ===== Move M/U to correct device/dtype =====
    if M_tensor is not None:
        M_tensor = M_tensor.to(device=device, dtype=dtype)
    if U_tensor is not None:
        U_tensor = U_tensor.to(device=device, dtype=dtype)

    # ===== Diffusion loop =====
    H = X_tensor
    if show_progress:
        progress_context = Progress(
            SpinnerColumn(),
            TextColumn("Diffusing"),
            BarColumn(),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            transient=True,
        )
        progress = progress_context.__enter__()
        task_id = progress.add_task("diffuse", total=n_steps)
    else:
        progress = None
        task_id = None

    for _ in range(n_steps):
        agg = torch.zeros_like(H)
        agg.index_add_(0, dst, H[src] * w.unsqueeze(-1))
        agg = agg / deg.unsqueeze(-1)

        # Feature coupling
        if M_tensor is not None:
            agg = agg @ M_tensor
        elif U_tensor is not None and beta != 0.0:
            # agg @ (I + beta U U^T) = agg + beta (agg U) U^T
            agg = agg + beta * ((agg @ U_tensor) @ U_tensor.T)

        H = (1 - alpha) * H + alpha * agg

        # Clamp anchored nodes
        if anchor_mask_tensor is not None:
            H[anchor_mask_tensor] = X_tensor[anchor_mask_tensor]

        if progress is not None:
            progress.advance(task_id)

    if progress is not None:
        progress_context.__exit__(None, None, None)

    # ===== Return =====
    if return_numpy:
        return H.cpu().numpy()
    else:
        return H


if __name__ == "__main__":
    # Simple test case
    X = torch.tensor([[1.0], [2.0], [3.0], [4.0]], dtype=torch.float32)
    edge_index = torch.tensor(
        [[0, 1, 2, 3, 0, 1], [1, 2, 3, 0, 2, 3]],
        dtype=torch.long,
    )
    X_smooth = diffuse(X, edge_index, alpha=0.5, n_steps=1)
    print("Original X:\n", X)
    print("Smoothed X:\n", X_smooth)
