library(tidyverse)
library(ggbrandon)
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

theme_set(theme_brandon(base_line_size = 0.3, base_size = 7))
ht_opt(
  COLUMN_ANNO_PADDING = unit(0.1, "mm"),
  ROW_ANNO_PADDING = unit(0.1, "mm"),
  TITLE_PADDING = unit(c(0.1, 0.1), "mm")
)
ht_opt$annotation_use_raster <- TRUE

# ----------------------------------------------------------------------
output_dir <- "figures/f4"
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
  annotation_label = c(""),
  show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = 5),
  simple_anno_size = unit(0.6, "lines"),
  raster_quality = 5
)

col_fun <- colorRamp2(
  seq(-2.5, 2.5, length.out = 8),
  rev(RColorBrewer::brewer.pal(8, "RdYlBu"))
)
# ----------------------------------------------------------------------
do_heatmap <- function(res, branch_name, show_legend = TRUE) {
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
    show_legend = FALSE,
    annotation_label = "Peak\nactivity",
    simple_anno_size = unit(0.4, "lines"),
    raster_quality = 5
  )

  Heatmap(
    mat,
    use_raster = TRUE,
    raster_quality = 5,
    col = col_fun,
    top_annotation = top_anno,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    row_order = row_order,
    show_heatmap_legend = FALSE,
    row_title = sprintf(
      "%s TDA TRs",
      scales::comma(nrow(mat))
    ),
    row_title_gp = gpar(fontsize = 6),
    column_title = branches_names[branch_name],
    column_title_gp = gpar(fontsize = 6),
    left_annotation = left_anno,
  )
}

# ----------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df()
res <- "box/results/lamian/z/tde.rds" %>%
  readRDS()

ht_padding <- unit(c(0, 0, 0, 0), "mm")

hmap1 <- grid.grabExpr(draw(
  do_heatmap(res, branches[[1]]),
  background = "transparent",
  padding = ht_padding
))

hmap2 <- grid.grabExpr(draw(
  do_heatmap(res, branches[[2]]),
  background = "transparent",
  padding = ht_padding
))

hmap3 <- grid.grabExpr(draw(
  do_heatmap(res, branches[[3]]),
  background = "transparent",
  padding = ht_padding
))

lgd_scaled <- Legend(
  col_fun = col_fun,
  title = "Scaled\nactivity",
  grid_width = unit(0.4, "lines"),
  legend_height = unit(1.5, "lines"),
  at = c(-2, 0, 2),
  labels = paste0("  ", c(-2, 0, 2)),
  title_gp = gpar(fontsize = 6),
  labels_gp = gpar(fontsize = 5)
)

lgd_peak <- Legend(
  labels = c("Early", "Middle", "Late"),
  legend_gp = gpar(fill = RColorBrewer::brewer.pal(3, "YlOrRd")),
  title = "Peak\nactivity",
  title_gp = gpar(fontsize = 6),
  labels_gp = gpar(fontsize = 5),
  grid_height = unit(0.4, "lines"),
  grid_width = unit(0.4, "lines")
)

lgd_pseudotime <- Legend(
  col_fun = colorRamp2(
    seq(0, 1000, length.out = 8),
    RColorBrewer::brewer.pal(8, "YlGnBu")
  ),
  title = "Pseudo-\ntime",
  grid_width = unit(0.4, "lines"),
  legend_height = unit(1.5, "lines"),
  at = c(1000, 500, 0),
  labels = paste0("  ", c("Late", "Middle", "Early")),
  title_gp = gpar(fontsize = 6),
  labels_gp = gpar(fontsize = 5)
)

lgd_all <- grid.grabExpr(draw(packLegend(
  lgd_scaled, lgd_pseudotime, lgd_peak,
  direction = "horizontal", gap = unit(2, "mm")
)))

combined <- wrap_elements(hmap1) +
  plot_spacer() +
  wrap_elements(hmap2) +
  plot_spacer() +
  wrap_elements(hmap3) +
  wrap_elements(lgd_all) +
  plot_layout(nrow = 1, widths = c(1, 0.07, 1, 0.07, 1, 0.7)) &
  theme(plot.margin = margin(b = 2))


save_figure(output_dir, "heatmaps_tde", plot = combined, width = 7.1, height = 1)
