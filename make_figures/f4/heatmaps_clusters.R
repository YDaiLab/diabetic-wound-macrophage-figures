library(tidyverse)
library(ggbrandon)
library(patchwork)
library(ggh4x)
library(legendry)
library(RColorBrewer)
library(Seurat)
library(ComplexHeatmap)
library(circlize)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 7))
ht_opt$annotation_use_raster <- TRUE

# ----------------------------------------------------------------------
output_dir <- "figures/f4"
source("make_figures/cfg.R")
source("make_figures/load_seurat.R")
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df()
res <- "results/bitfam_runs/pvals_combined.csv" %>%
  read_csv() %>%
  mutate(cluster = factor(cluster, levels(df$clusterid)))

genes <- res %>%
  arrange(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(order_by = avg_log2FC, by = cluster, n = 5) %>%
  pull(gene) %>%
  unique()

# ----------------------------------------------------------------------
cells <- load_cells()
cells <- ScaleData(cells, assay = "Z_avg", features = genes)
mat1 <- GetAssayData(cells, assay = "Z_avg", layer = "scale.data")

cells <- ScaleData(cells, assay = "RNA", features = genes)
mat2 <- GetAssayData(cells, assay = "RNA", layer = "scale.data")

all(colnames(mat1) == colnames(mat2))
# ----------------------------------------------------------------------

anno_df <- df %>%
  column_to_rownames("Row.names") %>%
  select(clusterid, sample)

anno_df <- anno_df[colnames(mat1), ]

top_anno <- HeatmapAnnotation(
  df = anno_df,
  col = list(
    clusterid = cluster_pal,
    sample = sample_pal
  ),
  gap = unit(0, "lines"),
  simple_anno_size = unit(0.3, "lines"),
  show_legend = c(FALSE, TRUE),
  annotation_legend_param = list(
    sample = list(
      title = "Sample",
      title_gp = gpar(fontsize = 6),
      labels_gp = gpar(fontsize = 6)
    )
  ),
  show_annotation_name = FALSE,
  raster_quality = 5
)

column_order <- df %>%
  arrange(clusterid, sample) %>%
  pull(Row.names)

col1 <- colorRamp2(
  breaks = seq(-2.5, 2.5, length.out = 8),
  colors = rev(brewer.pal(8, "RdYlBu"))
)

col2 <- colorRamp2(
  breaks = seq(-2.5, 2.5, length.out = 8),
  colors = rev(brewer.pal(8, "RdBu"))
)

ht_opt(
  COLUMN_ANNO_PADDING = unit(0.1, "mm"),
  TITLE_PADDING = unit(c(0.1, 0.1), "mm")
)

ht1 <- Heatmap(
  mat1,
  use_raster = TRUE,
  raster_quality = 5,
  col = col1,
  top_annotation = top_anno,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  column_order = column_order,
  row_order = genes,
  column_split = anno_df$clusterid,
  column_title_gp = gpar(fontsize = 6),
  column_gap = unit(0.1, "lines"),
  row_names_side = "left",
  heatmap_legend_param = list(
    title = "Scaled activity",
    grid_height = unit(0.5, "lines"),
    legend_width = unit(2, "lines"),
    direction = "horizontal",
    at = c(-2, 0, 2),
    labels = paste0("\n", c(-2, 0, 2)),
    title_gp = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 5)
  )
)

ht2 <- Heatmap(
  mat2[rownames(mat1), ],
  use_raster = TRUE,
  raster_quality = 5,
  col = col2,
  top_annotation = top_anno,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 5, fontface = "italic"),
  column_order = column_order,
  row_order = genes,
  column_split = anno_df$clusterid,
  column_title_gp = gpar(fontsize = 6),
  column_gap = unit(0.1, "lines"),
  row_names_side = "left",
  heatmap_legend_param = list(
    title = "Scaled expression",
    grid_height = unit(0.5, "lines"),
    legend_width = unit(2, "lines"),
    direction = "horizontal",
    at = c(-2, 0, 2),
    labels = paste0("\n", c(-2, 0, 2)),
    title_gp = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 5)
  )
)

# ----------------------------------------------------------------------
# Save each heatmap separately, with NO legends at all
# (no heatmap legend, no annotation legend)

hmap1 <- grid.grabExpr(
  draw(ht1,
    background = "transparent",
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE
  )
)

hmap2 <- grid.grabExpr(
  draw(ht2,
    background = "transparent",
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE
  )
)

# ----------------------------------------------------------------------
# Build a combined legend panel: annotation legend + both heatmap legends

lgd_sample <- Legend(
  labels = names(sample_pal),
  title = "Sample",
  legend_gp = gpar(fill = sample_pal),
  title_gp = gpar(fontsize = 6),
  labels_gp = gpar(fontsize = 5.5),
  grid_height = unit(0.5, "lines"),
  grid_width = unit(0.5, "lines")
)

lgd_activity <- Legend(
  col_fun = col1,
  title = "Scaled activity",
  at = c(-2, 0, 2),
  labels = paste0("\n", c(-2, 0, 2)),
  direction = "horizontal",
  grid_height = unit(0.4, "lines"),
  legend_width = unit(2, "lines"),
  title_gp = gpar(fontsize = 6),
  labels_gp = gpar(fontsize = 5)
)

lgd_expression <- Legend(
  col_fun = col2,
  title = "Scaled expression",
  at = c(-2, 0, 2),
  labels = paste0("\n", c(-2, 0, 2)),
  direction = "horizontal",
  grid_height = unit(0.4, "lines"),
  legend_width = unit(2, "lines"),
  title_gp = gpar(fontsize = 6),
  labels_gp = gpar(fontsize = 5)
)

leg <- grid.grabExpr(
  draw(packLegend(
    lgd_sample,
    lgd_activity,
    lgd_expression,
    direction = "vertical",
    gap = unit(1, "lines")
  ))
)

# ----------------------------------------------------------------------
# Save outputs

save_figure(output_dir, "heatmaps_clusters_hmap1", plot = hmap1, width = 3, height = 2.8)
save_figure(output_dir, "heatmaps_clusters_hmap2", plot = hmap2, width = 3, height = 2.8)
save_figure(output_dir, "heatmaps_clusters_legend", plot = leg, width = 1.08, height = 2.15)
