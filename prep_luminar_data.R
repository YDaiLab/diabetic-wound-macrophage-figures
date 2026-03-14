library(tidyverse)

# --- Configuration ---
factors <- c("Irf4", "Cebpa", "Mafb", "Kdm1a")
scenario <- "KO" # "KO" or "OE"

# --- Paths ---
scenario_dirs <- c(KO = "pgd_celloracle", OE = "pgd_celloracle_overexpress")
input_path <- file.path("box/results", scenario_dirs[[scenario]], "perturb_scores.parquet")
output_dir <- file.path("luminar_data", scenario)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- arrow::read_parquet(input_path) |>
  filter(is.na(branch))

# Pseudotime (same for any TF within a scenario)
df |>
  filter(factor == factors[[1]]) |>
  select(x, y, dx = ref_dx, dy = ref_dy) |>
  write_csv(file.path(output_dir, "pseudotime.csv"))

# Per-TF perturbation vectors
for (factor_ in factors) {
  df |>
    filter(factor == factor_) |>
    select(x, y, dx, dy) |>
    write_csv(file.path(output_dir, paste0(factor_, ".csv")))
}
