library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(ggarrow)
library(ggarchery)
library(ggrastr)

theme_set(theme_nature(base_size = fs_small))

# ----------------------------------------------------------------------
output_dir <- "figures/f1"
source("make_figures/cfg.R")
source("make_figures/load_edge_tbl.R")

axis <- legendry::guide_axis_base(cap = I(c(-Inf, 0.15)))
# ----------------------------------------------------------------------

df <- "results/pgd_iter_meta.csv" |>
  read_csv() |>
  mutate(clusterid = factor(clusterid))

tree <-
  "box/results/lamian/evaluate_uncertainty.rds" %>%
  readRDS()

edge_tbl <- df %>%
  group_split(alpha) %>%
  map(~ load_edge_tbl(tree, .x) %>% mutate(alpha = .x$alpha[1])) %>%
  bind_rows() %>%
  arrange(alpha, detection.rate)

p <- df |>
  filter(alpha %in% c(0, 0.6)) |>
  ggplot(aes(umap_1, umap_2)) +
  geom_point_rast(aes(color = pseudotime), size = 0.35, stroke = 0, raster.dpi = 600) +
  geom_arrow_segment(
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    edge_tbl %>%
      filter(detection.rate > 0.40, alpha %in% c(0, 0.6)),
    lwd = 0.3,
    position = position_attractsegment(
      start_shave = 0.1,
      end_shave = 0.1
    )
  ) +
  facet_wrap(~alpha,
    ncol = 1, axes = "all", scales = "free",
    labeller = label_bquote(alpha == .(alpha))
  ) +
  scale_color_distiller(
    palette = "YlGnBu",
    direction = 1,
    guide = guide_colorbar(
      theme = theme(
        legend.key.width = unit(0.3, "lines"),
        legend.key.height = unit(1, "lines"),
        legend.text = element_text(size = fs_min),
        legend.ticks = element_blank(),
        legend.title = element_text(size = fs_min, hjust = 0.5, margin = margin(b = 7)),
        legend.title.position = "top",
      )
    ),
    breaks = range(df$pseudotime),
    labels = c("Early", "Late")
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Pseudotime"
  ) +
  theme_embedding(arrow = arrow, axis = axis) +
  theme(
    axis.line = element_blank(),
    panel.spacing.y = unit(2, "lines"),
    legend.position = "right",
    legend.justification = "top",
    legend.direction = "vertical",
    strip.text = element_blank(),
    plot.margin = margin(),
    legend.margin = margin()
  )

save_figure(output_dir, "umaps_pgd_iter", plot = p, width = 1.35, height = 2.5)
