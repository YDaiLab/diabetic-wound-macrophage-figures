library(tidyverse)
library(Seurat)

# ----------------------------------------------------------------------
reference_file <- "box/results/seurat/cells.rds"
run0_file <- "box/results/bitfam/results/acts.rds"
run1_file <- "box/results/bitfam_runs/run_1/acts.rds"
run2_file <- "box/results/bitfam_runs/run_2/acts.rds"
run3_file <- "box/results/bitfam_runs/run_3/acts.rds"

output_file <- "results/bitfam_runs/pvals.csv"
# ----------------------------------------------------------------------

reference <-
  reference_file %>%
  readRDS()

extract_pvals_from_cells <- function(cells) {
  require(tibble)
  require(dplyr)

  cells$clusterid <- reference$clusterid[colnames(cells)]
  Idents(cells) <- "clusterid"
  DefaultAssay(cells) <- "Z"

  FindAllMarkers(
    cells,
    min.pct = 0,
    logfc.threshold = 0,
    return.thresh = Inf,
    only.pos = FALSE
  ) %>%
    as_tibble()
}

pvals <-
  run0_file %>%
  readRDS() %>%
  extract_pvals_from_cells()

# Runs
pvals1 <-
  run1_file %>%
  readRDS() %>%
  extract_pvals_from_cells()

pvals2 <-
  run2_file %>%
  readRDS() %>%
  extract_pvals_from_cells()

pvals3 <-
  run3_file %>%
  readRDS() %>%
  extract_pvals_from_cells()

df <- bind_rows(
  pvals %>%
    mutate(run = "Original"),
  pvals1 %>%
    mutate(run = "Run 1"),
  pvals2 %>%
    mutate(run = "Run 2"),
  pvals3 %>%
    mutate(run = "Run 3")
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(df, output_file)
