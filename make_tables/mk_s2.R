source("make_tables/cfg.R")

file.copy(
  "results/enrichment_results/GSEA.csv",
  get_table_path(2, "cluster_gsea_enrichment.csv"),
  overwrite = TRUE
)
