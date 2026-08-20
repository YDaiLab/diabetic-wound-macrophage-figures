library(tidyverse)

table_base_path <- file.path(
  "/Users/brandonlukas/Library/CloudStorage/Box-Box",
  "Brandon Lukas",
  "Koh",
  "Manuscript and Figures",
  "NPJ files",
  "Supplementary Information"
)

# kind = "Data" for spreadsheet-scale files, "Table" for page-formatted ones
get_table_path <- function(table_num, filename, kind = "Data") {
  name <- tools::file_path_sans_ext(filename)
  ext <- tools::file_ext(filename)
  out <- sprintf("Supplementary_%s_%d_%s.%s", kind, table_num, name, ext)
  path <- file.path(table_base_path, out)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
