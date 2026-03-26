source("make_tables/cfg.R")

file.copy(
  "results/bitfam_runs/pvals_combined.csv",
  get_table_path("Table S4", "pvals_combined.csv"),
  overwrite = TRUE
)
