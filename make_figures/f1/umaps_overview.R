library(tidyverse)
library(ggbrandon)
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(ggarrow)
library(ggarchery)
library(ggrastr)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 6))

# ----------------------------------------------------------------------
output_dir <- "figures/f1"
source("make_figures/cfg.R")
source("make_figures/load_edge_tbl.R")

axis <- legendry::guide_axis_base(cap = I(c(-Inf, 0.15)))
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df()
df2 <- load_df(pgd = TRUE)

tree <-
  "box/results/lamian/evaluate_uncertainty.rds" %>%
  readRDS()

edge_tbl <- load_edge_tbl(tree, df) %>%
  arrange(detection.rate)

edge_tbl2 <- load_edge_tbl(tree, df2) %>%
  arrange(detection.rate)

p <- df %>%
  ggplot(aes(umap_1, umap_2)) +
  geom_point_rast(aes(color = sample), size = 0.35, stroke = 0, raster.dpi = 600) +
  scale_color_manual(values = sample_pal) +
  labs(
    x = NULL,
    y = NULL,
    color = NULL
  ) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.2),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  theme_dimred2(arrow = arrow, axis = axis) +
  theme(axis.line = element_blank())

save_figure(output_dir, "umaps_sample", plot = p, width = 0.9, height = 0.9)

p <- df %>%
  ggplot(aes(umap_1, umap_2)) +
  geom_point_rast(aes(color = clusterid), size = 0.35, stroke = 0, raster.dpi = 600) +
  geom_arrow_segment(
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    edge_tbl %>%
      filter(detection.rate > 0.40),
    lwd = 0.3,
    position = position_attractsegment(
      start_shave = 0.1,
      end_shave = 0.1
    )
  ) +
  scale_color_manual(values = cluster_pal) +
  labs(
    x = NULL,
    y = NULL,
    color = NULL
  ) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.2),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  theme_dimred2(arrow = arrow, axis = axis) +
  theme(axis.line = element_blank())

save_figure(output_dir, "umaps_clusters_trajectory", plot = p, width = 0.9, height = 0.9)
