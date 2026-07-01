library(tidyverse)

# --- Configuration ---
factors <- c("Irf4", "Cebpa", "Mafb", "Kdm1a")
scenario <- "KO" # "KO" or "OE"

# --- Paths ---
# Two CellOracle layouts: PGD embedding and the standard UMAP embedding.
pgd_dirs <- c(KO = "pgd_celloracle", OE = "pgd_celloracle_overexpress")
umap_dirs <- c(KO = "celloracle", OE = "celloracle_overexpress")
pgd_path <- file.path("box/results", pgd_dirs[[scenario]], "perturb_scores.parquet")
umap_path <- file.path("box/results", umap_dirs[[scenario]], "perturb_scores.parquet")
output_dir <- file.path("luminar_data", scenario)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- arrow::read_parquet(pgd_path) |>
  filter(is.na(branch))

# Pseudotime (reference development flow) on the PGD layout.
# ref_dx/ref_dy are the same for any TF within a scenario.
df |>
  filter(factor == factors[[1]]) |>
  select(x, y, dx = ref_dx, dy = ref_dy) |>
  write_csv(file.path(output_dir, "pseudotime.csv"))

# Pseudotime (reference development flow) on the standard UMAP layout.
arrow::read_parquet(umap_path) |>
  filter(is.na(branch)) |>
  filter(factor == factors[[1]]) |>
  select(x, y, dx = ref_dx, dy = ref_dy) |>
  write_csv(file.path(output_dir, "umap_pseudotime.csv"))

# Per-TF perturbation vectors (PGD layout)
for (factor_ in factors) {
  df |>
    filter(factor == factor_) |>
    select(x, y, dx, dy) |>
    write_csv(file.path(output_dir, paste0(factor_, ".csv")))
}
