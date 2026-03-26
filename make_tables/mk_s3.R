source("make_tables/cfg.R")
library(Lamian)

df1 <- read_csv("results/lamian/rna_tde.csv")
df2 <- read_csv("results/lamian/rna_tde_go.csv")

sheets <- list(
  "TDE" = df1,
  "TDE_Pathway_Enrichment" = df2
)

output_file <- get_table_path("Table S3", "lamian_tde.xlsx")
writexl::write_xlsx(sheets, output_file)
