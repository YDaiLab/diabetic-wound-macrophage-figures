source("make_tables/cfg.R")

file.copy(
  "results/de_tests.xlsx",
  get_table_path("Table S1", "de_tests.xlsx"),
  overwrite = TRUE
)
