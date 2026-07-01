source("make_tables/cfg.R")

file.copy(
  "results/bitfam_runs/pvals_combined.csv",
  get_table_path(4, "cluster_differential_tr_activity.csv"),
  overwrite = TRUE
)
