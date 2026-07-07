library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(ggrastr)

theme_set(theme_nature())

# ----------------------------------------------------------------------
output_dir <- "figures/f2"
source("make_figures/cfg.R")
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df()

p1 <- df %>%
  ggplot(aes(umap_1, umap_2)) +
  geom_point_rast(aes(color = sample), size = 0.3, stroke = 0, raster.dpi = 600) +
  scale_color_manual(values = sample_pal) +
  theme_embedding(arrow = arrow) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    color = NULL
  ) +
  guides(
    color = guide_legend_cross(
      key = key_group_split(sep = " "),
      override.aes = list(size = 3, shape = 15),
      swap = TRUE,
      col_text = element_text(margin = margin(b = 2))
    )
  ) +
  theme(
    axis.title = element_text(size = fs_small),
    legend.text.position = c("top", "left"),
    legend.position = "top",
    legend.justification = "left",
    legend.location = "plot",
    legend.margin = margin(),
    legend.key.spacing = rel(0),
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(5, "mm")
  )

p2 <- df %>%
  ggplot(aes(umap_1, umap_2)) +
  geom_point_rast(aes(color = clusterid), size = 0.3, stroke = 0, raster.dpi = 600) +
  geom_shadowtext(
    aes(label = clusterid),
    df %>%
      group_by(clusterid) %>%
      summarise_if(is.numeric, median),
    color = "#000",
    bg.color = "#fff",
    size = pt2mm(fs_base)
  ) +
  scale_color_manual(values = cluster_pal) +
  theme_embedding(arrow = arrow) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    color = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5)
    )
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = fs_small)
  )

combined <- p1 + plot_spacer() + p2 +
  plot_layout(nrow = 1, widths = c(1, 0.05, 1)) & theme(plot.margin = margin())

save_figure(output_dir, "umaps", plot = combined, width = 3.7, height = 2.1)
