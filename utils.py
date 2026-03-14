import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import default_converter, numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.rlike.container import NamedList


def load_orders_from_lamian(path: str) -> dict:
    # Load Lamian infer_tree result
    with localconverter(default_converter + numpy2ri.converter + pandas2ri.converter):
        result = ro.r(f"readRDS('{path}')")

    orders = _named_list_to_dict(_named_list_to_dict(result)["order"])
    return orders


def run_umap_seurat(embeddings: np.ndarray) -> pd.DataFrame:
    string = """
local({
  library(tidyverse)
  library(Seurat)

  X <- as.matrix(X)

  cells <- CreateSeuratObject(counts = t(X))
  cells@reductions[["pca"]] <- CreateDimReducObject(
    embeddings = X,
    assay = "RNA"
  )
  cells <- cells %>%
    FindNeighbors() %>%
    RunUMAP(dims = 1:30)

  cells@reductions[["umap"]]@cell.embeddings
})
"""
    with localconverter(default_converter + numpy2ri.converter + pandas2ri.converter):
        ro.globalenv["X"] = embeddings
        umap_result = ro.r(string)
        ro.globalenv.clear()

    df = pd.DataFrame(umap_result, index=embeddings.index, columns=["UMAP_1", "UMAP_2"])
    return df


def _named_list_to_dict(named_list: NamedList) -> dict:
    return dict(
        zip(
            list(map(str, named_list._NamedList__names)),
            list(named_list._NamedList__list),
        )
    )
