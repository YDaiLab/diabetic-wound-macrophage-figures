library(tidyverse)

table_base_path <- file.path(
  "/Users/brandonlukas/Library/CloudStorage/Box-Box",
  "Brandon Lukas",
  "Koh",
  "Manuscript and Figures",
  "Supplementary Data"
)

get_table_path <- function(table_num, filename) {
  name <- tools::file_path_sans_ext(filename)
  ext <- tools::file_ext(filename)
  out <- sprintf("Supplementary_Data_%d_%s.%s", table_num, name, ext)
  path <- file.path(table_base_path, out)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
