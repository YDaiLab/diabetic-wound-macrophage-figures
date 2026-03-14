library(tidyverse)
library(Seurat)

# ----------------------------------------------------------------------
run0_file <- "box/results/bitfam/results/acts.rds"
run1_file <- "box/results/bitfam_runs/run_1/acts.rds"
run2_file <- "box/results/bitfam_runs/run_2/acts.rds"
run3_file <- "box/results/bitfam_runs/run_3/acts.rds"

output_file <- "results/bitfam_runs/acts.csv"
# ----------------------------------------------------------------------

extract_acts_from_cells <- function(cells) {
  require(tibble)
  require(tidyr)
  require(dplyr)

  data.frame(
    cells@assays[["Z"]]@data,
    row.names = rownames(cells@assays[["Z"]]@data)
  ) %>%
    rownames_to_column("gene") %>%
    pivot_longer(where(is.numeric))
}

acts <-
  run0_file %>%
  readRDS() %>%
  extract_acts_from_cells()

# Runs
acts1 <-
  run1_file %>%
  readRDS() %>%
  extract_acts_from_cells()

acts2 <-
  run2_file %>%
  readRDS() %>%
  extract_acts_from_cells()

acts3 <-
  run3_file %>%
  readRDS() %>%
  extract_acts_from_cells()

df <- bind_rows(
  acts %>%
    mutate(run = "Original"),
  acts1 %>%
    mutate(run = "Run 1"),
  acts2 %>%
    mutate(run = "Run 2"),
  acts3 %>%
    mutate(run = "Run 3")
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(df, output_file)
