library(tidyverse)
library(Seurat)

load_cells <- function() {
  readRDS("box/results/seurat/cells.rds")
}
