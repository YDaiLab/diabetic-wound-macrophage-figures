source("make_tables/cfg.R")

file.copy(
  "results/de_tests.xlsx",
  get_table_path(1, "cluster_differential_expression.xlsx"),
  overwrite = TRUE
)
