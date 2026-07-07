library(tidyverse)
source("assets/theme_nature.R")
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(ggnewscale)
library(ggarrow)
library(ggarchery)
library(Lamian)
library(ComplexHeatmap)
library(circlize)

theme_set(theme_nature())
ht_opt$annotation_use_raster <- TRUE

# ----------------------------------------------------------------------
output_dir <- "figures/f3"
source("make_figures/cfg.R")
source("make_figures/load_edge_tbl.R")

branches <- c(
  "branch: 12,6,7,5",
  "branch: 12,6,9,3,10",
  "branch: 12,6,9,3,4,8"
)

branches_names <- c(
  "branch 1: 12,6,7,5",
  "branch 2: 12,6,9,3,10",
  "branch 3: 12,6,9,3,4,8"
)
branches_names <- setNames(branches_names, branches)
# ----------------------------------------------------------------------
top_anno <- HeatmapAnnotation(
  pseudotime = 1:1000,
  col = list(
    pseudotime = colorRamp2(
      seq(0, 1000, length.out = 8),
      RColorBrewer::brewer.pal(8, "YlGnBu")
    )
  ),
  annotation_label = c("Pseudotime"),
  show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = fs_min),
  simple_anno_size = unit(0.8, "lines"),
  raster_quality = 5
)

col_fun <- colorRamp2(
  seq(-2, 2, length.out = 8),
  rev(RColorBrewer::brewer.pal(8, "RdBu"))
)
# ----------------------------------------------------------------------
do_heatmap <- function(res, branch_name) {
  require(grid)
  require(dplyr)

  res <- res[[branch_name]]
  diffgene <- rownames(res$statistics)[res$statistics[, 1] < 0.05]
  tmp <-
    getPopulationFit(res, gene = diffgene, type = "time")

  mat <- t(scale(t(tmp)))
  peak_time <- apply(mat, 1, which.max)
  row_order <- names(sort(peak_time))

  left_anno <- rowAnnotation(
    cluster = factor(
      case_when(
        peak_time < 10 ~ "Early",
        peak_time > 990 ~ "Late",
        .default = "Middle"
      ),
      c("Early", "Middle", "Late")
    ),
    col = list(
      cluster = setNames(
        RColorBrewer::brewer.pal(3, "YlOrRd"),
        c("Early", "Middle", "Late")
      )
    ),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      title_gp = gpar(fontsize = fs_small),
      labels_gp = gpar(fontsize = fs_small),
      grid_height = unit(0.4, "lines"),
      grid_width = unit(0.4, "lines")
    ),
    annotation_label = "Peak\nexpression",
    simple_anno_size = unit(0.6, "lines"),
    raster_quality = 5
  )

  Heatmap(
    mat,
    use_raster = TRUE,
    raster_quality = 5,
    left_annotation = left_anno,
    col = col_fun,
    top_annotation = top_anno,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    row_order = row_order,
    heatmap_legend_param = list(
      title = "Scaled\nexpression",
      grid_width = unit(0.5, "lines"),
      legend_height = unit(1.5, "lines"),
      at = c(-2, 0, 2),
      labels = paste0("  ", c(-2, 0, 2)),
      title_gp = gpar(fontsize = fs_small),
      labels_gp = gpar(fontsize = fs_small)
    ),
    row_title = sprintf(
      "%s TDE genes",
      scales::comma(nrow(mat)),
      branch_name
    ),
    row_title_gp = gpar(fontsize = fs_small),
    column_title = branches_names[branch_name],
    column_title_gp = gpar(fontsize = fs_base)
  )
}

do_heatmap_save <- function(res, branch_name, show_legend = TRUE) {
  require(ggplot2)

  ht <- do_heatmap(res, branch_name)
  hmap <-
    grid.grabExpr(draw(ht,
      background = "transparent",
      show_heatmap_legend = show_legend,
      show_annotation_legend = show_legend
    ))

  save_figure(output_dir, sprintf("heatmaps_tde_%s", branch_name),
    plot = hmap,
    width = 3, height = 1.8
  )
}
# ----------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df()
res <- "box/results/lamian/rna/tde.rds" %>%
  readRDS()

do_heatmap_save(res, branches[[1]])
do_heatmap_save(res, branches[[2]])
do_heatmap_save(res, branches[[3]])
