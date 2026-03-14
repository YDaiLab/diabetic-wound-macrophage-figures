library(tidyverse)
library(ggbrandon)
library(patchwork)
library(ggh4x)
library(legendry)
library(RColorBrewer)
library(Seurat)
library(ComplexHeatmap)
library(circlize)

ht_opt$annotation_use_raster <- TRUE

theme_set(theme_brandon(base_line_size = 0.3, base_size = 7))

# ----------------------------------------------------------------------
output_dir <- "figures/f1"
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
  simple_anno_size = unit(0.5, "lines"),
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

ht1 <- Heatmap(
  mat1,
  use_raster = TRUE,
  raster_quality = 5,
  col = col1,
  top_annotation = top_anno,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  column_order = column_order,
  row_order = genes,
  column_split = anno_df$clusterid,
  column_title = NULL,
  row_title = "Top TRs per cluster",
  row_title_gp = gpar(fontsize = 5),
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

hmap <- grid.grabExpr(
  draw(ht1,
    background = "transparent",
    merge_legend = TRUE,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE
  )
)

save_figure(output_dir, "heatmap_tf_activity", plot = hmap, width = 2.2, height = 1.9)
