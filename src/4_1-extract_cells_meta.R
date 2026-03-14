library(tidyverse)
library(Seurat)

# ----------------------------------------------------------------------
input_file <- "box/results/seurat/cells.rds"
output_file <- "results/cells_meta.csv"
# ----------------------------------------------------------------------


cells <- readRDS(input_file)

df <- merge(
  cells@reductions$umap@cell.embeddings,
  cells@meta.data,
  by = 0
) %>%
  as_tibble() %>%
  mutate(
    sample = case_match(
      orig.ident,
      "WT_D3" ~ "Non-diabetic D3",
      "WT_D6" ~ "Non-diabetic D6",
      "WT_D10" ~ "Non-diabetic D10",
      "db_D3" ~ "Diabetic D3",
      "db_D6" ~ "Diabetic D6",
      "db_D10" ~ "Diabetic D10"
    ),
    sample = fct(sample, c(
      "Non-diabetic D3",
      "Non-diabetic D6",
      "Non-diabetic D10",
      "Diabetic D3",
      "Diabetic D6",
      "Diabetic D10"
    ))
  ) %>%
  separate(sample,
    into = c("condition", "timepoint"),
    sep = " ",
    remove = FALSE
  ) %>%
  mutate(
    condition = fct(condition, c("Non-diabetic", "Diabetic")),
    timepoint = fct(timepoint, c("D3", "D6", "D10"))
  )

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(df, output_file)
