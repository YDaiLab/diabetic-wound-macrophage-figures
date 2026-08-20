library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
source("assets/theme_nature.R")

# Revision R2-2 — BIRDccNEST alternative-method cross-check (panel j). The
# companion to the cell graph (panel h) and oMST trajectory (panel i): how the
# BIRDccNEST communities map onto the 12 Lamian clusters. Row-normalized overlap
# (fraction of each community's cells falling in each Lamian cluster).
#
# ComplexHeatmap: rows (communities) carry a color strip in the SAME community
# palette as panels h/i; columns (Lamian clusters) carry a color strip in the
# SAME cluster_pal as panels b-f/h/i-lamian. The colored strips let the reader
# tie each row/column to the other panels by color instead of chasing numbers.
# Communities (rows) ordered by mean pseudotime (earliest on top); clusters
# (columns) by flow-weighted community pseudotime, so shared structure reads as a
# diagonal.
# ----------------------------------------------------------------------
output_dir <- "figures/revision/trajectory/supp_panels"
source("make_figures/cfg.R")

bdir <- "box/results/revision/birdccnest/joint_res1p5"
MIN_CELLS <- 100 # res1.5: keep 5 communities >=100 cells (drops 7=66, 1=50, 6=46)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# community palette (matches panels h/i)
community_pal <- c(
  "0" = "#4E79A7", "1" = "#F28E2B", "2" = "#E15759", "3" = "#76B7B2",
  "4" = "#59A14F", "5" = "#EDC948", "6" = "#B07AA1", "7" = "#FF9DA7",
  "8" = "#9C755F", "9" = "#BAB0AC", "10" = "#86BCB6", "11" = "#D37295",
  "12" = "#A0CBE8"
)

pt <- read_csv(file.path(bdir, "community_mean_pseudotime.csv"), show_col_types = FALSE)

# community x clusterid counts -> long, row-normalized; drop communities < MIN_CELLS
long <- read_csv(file.path(bdir, "community_overlap_clusterid.csv"), show_col_types = FALSE) %>%
  pivot_longer(-community, names_to = "clusterid", values_to = "n") %>%
  group_by(community) %>%
  filter(sum(n) >= MIN_CELLS) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  left_join(pt, by = "community")

# row order: communities by ascending pseudotime (earliest on top for ComplexHeatmap)
comm_levels <- pt %>%
  filter(community %in% unique(long$community)) %>%
  arrange(mean_pseudotime) %>%
  pull(community) %>%
  as.character()
# column order: clusters by flow-weighted community pseudotime (diagonal)
clus_levels <- long %>%
  group_by(clusterid) %>%
  summarise(o = weighted.mean(mean_pseudotime, n), .groups = "drop") %>%
  arrange(o) %>%
  pull(clusterid) %>%
  as.character()

# communities x clusters fraction matrix, ordered
mat <- long %>%
  transmute(community = as.character(community), clusterid = as.character(clusterid), frac) %>%
  pivot_wider(names_from = clusterid, values_from = frac) %>%
  column_to_rownames("community") %>%
  as.matrix()
mat <- mat[comm_levels, clus_levels, drop = FALSE]

# color strips tie rows/cols to the community (h/i) and cluster (b-f) palettes
left_anno <- rowAnnotation(
  community = rownames(mat), col = list(community = community_pal),
  simple_anno_size = unit(1.8, "mm"), show_legend = FALSE,
  show_annotation_name = FALSE
)
top_anno <- HeatmapAnnotation(
  `Lamian cluster` = colnames(mat), col = list(`Lamian cluster` = cluster_pal),
  simple_anno_size = unit(1.8, "mm"), show_legend = FALSE,
  show_annotation_name = TRUE, annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = fs_min)
)

col_fun <- colorRamp2(c(0, 1), c("#FFFFFF", "#0A574C"))

ht <- Heatmap(
  mat,
  name = "Fraction of\ncommunity", col = col_fun, na_col = "white",
  cluster_rows = FALSE, cluster_columns = FALSE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  border = TRUE, border_gp = gpar(col = "black", lwd = 0.6),
  left_annotation = left_anno, top_annotation = top_anno,
  row_names_side = "left", column_names_side = "top",
  row_names_gp = gpar(fontsize = fs_min), column_names_gp = gpar(fontsize = fs_min),
  column_names_rot = 0, column_names_centered = TRUE,
  row_title = "BIRDccNEST community", row_title_gp = gpar(fontsize = fs_min),
  width = ncol(mat) * unit(3.6, "mm"), height = nrow(mat) * unit(3.6, "mm"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- mat[i, j]
    if (!is.na(v) && v >= 0.15) {
      grid.text(sprintf("%.2f", v), x, y,
        gp = gpar(fontsize = fs_min - 1, col = if (v >= 0.5) "white" else "grey20")
      )
    }
  },
  heatmap_legend_param = list(
    at = c(0, 0.5, 1), title_gp = gpar(fontsize = fs_min, fontface = "plain"),
    labels_gp = gpar(fontsize = fs_min),
    grid_width = unit(2.5, "mm"), legend_height = unit(12, "mm")
  )
)

# Capture the ComplexHeatmap as a grid grob so it flows through the shared
# save_figure() (transparent background, no in-panel legends — the fraction scale
# and color strips are explained in the figure legend, as with panels b-f/h/i).
hmap <- grid.grabExpr(
  draw(ht,
    background = "transparent",
    merge_legend = TRUE,
    show_annotation_legend = FALSE
  )
)

# tiles fixed at 3.6 mm; canvas wide enough to contain the axis titles/labels
# (the grabbed grob draws at native size, so a too-small canvas clips it).
save_figure(output_dir, "panel_J_birdccnest_heatmap", plot = hmap, width = 3.1, height = 1.4)
