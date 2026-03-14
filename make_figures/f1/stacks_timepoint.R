library(tidyverse)
library(ggbrandon)
library(ComplexHeatmap)
library(circlize)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 6))

# ----------------------------------------------------------------------
output_dir <- "figures/f1"
source("make_figures/cfg.R")
# ----------------------------------------------------------------------

width <- 0.06
height <- 0.06 * 2
padding <- 0.01

do_timepoint_stack <- function(green, orange, suffix) {
  require(grid)
  require(ggplot2)

  mat <- matrix(c(1, 0), nrow = 2, ncol = 1)
  col <- colorRamp2(c(0, 1), c(orange, green))

  ht <- Heatmap(
    mat,
    col = col,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_heatmap_legend = FALSE,
    rect_gp = gpar(col = "black", lwd = 0.8),
    width = unit(width, "in"),
    height = unit(height, "in")
  )
  hmap <- grid.grabExpr(draw(ht, background = "transparent"))
  save_figure(output_dir, sprintf("stacks_timepoint_%s", suffix),
    plot = hmap,
    width = width + padding, height = height + padding
  )
}

do_timepoint_stack(sample_pal[["Non-diabetic D3"]], sample_pal[["Diabetic D3"]], "d3")
do_timepoint_stack(sample_pal[["Non-diabetic D6"]], sample_pal[["Diabetic D6"]], "d6")
do_timepoint_stack(sample_pal[["Non-diabetic D10"]], sample_pal[["Diabetic D10"]], "d10")
