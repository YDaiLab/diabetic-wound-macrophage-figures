source("make_tables/cfg.R")
library(Seurat)

factor_list <- "results/bitfam_runs/acts.csv" |>
  read_csv() |>
  pull(gene) |>
  unique()

cells <- readRDS("box/results/seurat/cells.rds")
dat <- AggregateExpression(
  cells,
  assays = "RNA",
  features = factor_list,
  return.seurat = TRUE
)

df <- GetAssayData(dat) |>
  as.data.frame() |>
  rownames_to_column("tf_name") |>
  pivot_longer(
    where(is.numeric),
    names_to = "cluster",
    values_to = "lognorm_pseudobulk_expression"
  ) |>
  mutate(cluster = gsub("g", "", cluster))

output_file <- get_table_path("Table S6", "pseudobulk_tf_expression.csv")
write_csv(df, output_file)
