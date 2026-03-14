library(tidyverse)
library(ggbrandon)
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(ggnewscale)
library(ggarrow)
library(ggarchery)
library(ggrastr)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 7))

# ----------------------------------------------------------------------
output_dir <- "figures/f3"
source("make_figures/cfg.R")
source("make_figures/load_edge_tbl.R")
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df(pgd = TRUE)

tree <-
  "box/results/lamian/evaluate_uncertainty.rds" %>%
  readRDS()

edge_tbl <- load_edge_tbl(tree, df) %>%
  arrange(detection.rate)

p1 <- df %>%
  ggplot(aes(umap_1, umap_2)) +
  geom_point_rast(
    aes(color = clusterid),
    size = 0.3,
    stroke = 0,
    show.legend = FALSE,
    raster.dpi = 600
  ) +
  scale_color_manual(values = cluster_pal) +
  new_scale_color() +
  geom_arrow_segment(
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = detection.rate
    ),
    edge_tbl,
    lwd = 0.5,
    position = position_attractsegment(
      start_shave = 0.2,
      end_shave = 0.2
    )
  ) +
  geom_shadowtext(
    aes(label = clusterid),
    df %>%
      group_by(clusterid) %>%
      summarise_if(is.numeric, median),
    color = "#000",
    bg.color = "#fff",
    size = 7 / .pt
  ) +
  colorspace::scale_color_continuous_sequential(
    palette = "Grays", rev = TRUE, begin = 0.3,
    breaks = scales::pretty_breaks(n = 3),
    labels = scales::percent_format(accuracy = 1),
    guide = guide_colorbar(
      theme = theme(
        legend.key.height = unit(1.5, "lines"),
        legend.key.width = unit(0.4, "lines"),
      )
    )
  ) +
  theme_dimred2(arrow = arrow) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    color = "Detection rate"
  ) +
  theme(
    axis.title = element_text(size = 6),
    legend.title = element_text(size = 6.75, hjust = 0.5, margin = margin(b = 7)),
    legend.text = element_text(size = 6),
    legend.position = "inside",
    legend.justification = c(-0.2, 1.4),
    legend.margin = margin(0, 0, 0, 0)
  )

p2 <- df %>%
  ggplot(aes(umap_1, umap_2)) +
  geom_point_rast(aes(color = pseudotime), size = 0.3, stroke = 0, raster.dpi = 600) +
  scale_color_distiller(
    palette = "YlGnBu",
    direction = 1,
    breaks = c(range(df$pseudotime), median(df$pseudotime)),
    labels = c("Early", "Late", "Middle"),
    guide = guide_colorbar(
      reverse = TRUE,
      theme = theme(
        legend.key.height = unit(1.5, "lines"),
        legend.key.width = unit(0.4, "lines")
      )
    )
  ) +
  labs(color = "Pseudotime") +
  new_scale_color() +
  geom_arrow_segment(
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = detection.rate
    ),
    edge_tbl,
    lwd = 0.5,
    position = position_attractsegment(
      start_shave = 0.2,
      end_shave = 0.2
    ),
    show.legend = FALSE
  ) +
  geom_shadowtext(
    aes(label = clusterid),
    df %>%
      group_by(clusterid) %>%
      summarise_if(is.numeric, median),
    color = "#000",
    bg.color = "#fff",
    size = 7 / .pt
  ) +
  colorspace::scale_color_continuous_sequential(
    palette = "Grays", rev = TRUE, begin = 0.3
  ) +
  theme_dimred2(arrow = arrow) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    color = "Detection rate"
  ) +
  theme(
    axis.title = element_text(size = 6),
    legend.title = element_text(size = 6.75, hjust = 0.5, margin = margin(b = 7)),
    legend.text = element_text(size = 6),
    legend.position = "inside",
    legend.justification = c(-0.2, 1.4),
    legend.margin = margin(0, 0, 0, 0)
  )

design <- "
444
123
"
combined <- p1 + plot_spacer() + p2 + plot_spacer() +
  plot_layout(
    design = design,
    widths = c(1, 0.10, 1),
    heights = c(0.3, 1)
  ) & theme(plot.margin = margin())

save_figure(output_dir, "umaps", plot = combined, width = 3.9, height = 2.1)
