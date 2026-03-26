library(tidyverse)

table_base_path <- file.path(
  "/Users/brandonlukas/Library/CloudStorage/Box-Box",
  "Brandon Lukas",
  "Koh",
  "Manuscript and Figures",
  "Tables"
)

get_table_path <- function(table_subdir, filename) {
  path <- file.path(table_base_path, table_subdir, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
