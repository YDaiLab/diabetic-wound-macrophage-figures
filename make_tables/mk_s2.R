source("make_tables/cfg.R")

file.copy(
  "results/enrichment_results/GSEA.csv",
  get_table_path("Table S2", "GSEA.csv"),
  overwrite = TRUE
)
