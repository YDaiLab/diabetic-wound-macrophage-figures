library(tidyverse)
library(ggbrandon)
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(ggnewscale)
library(ggarrow)
library(ggarchery)
library(colorspace)
library(ggrastr)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 5))

# ----------------------------------------------------------------------
output_dir <- "figures/f1"
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

mc_steps <- c(0, 10, 100, 500)
mc_dfs <- map(mc_steps, ~ {
  arrow::read_parquet(sprintf("box/results/pgd_celloracle/mc_transitions/%s.parquet", .x))
})
df2c <- arrow::read_parquet("box/results/pgd_celloracle/perturb_scores.parquet")

foi <- "Cebpa"
fois <- c("Cebpa")
dat2 <- map2_dfr(mc_dfs, mc_steps, ~ {
  .x %>%
    filter(factor == foi) %>%
    mutate(
      step = .y,
      factor = if (.y == 0) "Start" else sprintf("*%s* KO", foi)
    )
}) %>%
  mutate(
    step = sprintf("Step %s", step),
    step = fct_inorder(step),
    factor = fct_inorder(factor)
  )

p1 <- dat2 %>%
  ggplot(aes(umap_1, umap_2)) +
  geom_hex() +
  facet_nested(
    ~ factor + step,
    axes = "all",
    strip = strip_nested(
      background_x = list(
        element_rect(fill = "#EEEEEE"),
        element_rect(fill = "transparent")
      ),
      by_layer_x = TRUE
    )
  ) +
  scale_fill_distiller(
    palette = "Blues",
    direction = 1,
    guide = guide_colorbar(
      theme = theme(
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(0.3, "lines")
      )
    ),
    trans = scales::log10_trans(),
    labels = scales::comma_format(),
    breaks = c(1, 10, 100, 1000)
  ) +
  theme_dimred2(arrow = arrow) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    fill = "Counts"
  ) +
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = 5, vjust = 1),
    legend.position = "bottom",
    legend.title.position = "left",
    legend.justification = "left",
    legend.text.position = "bottom",
    legend.margin = margin(),
    strip.text = ggtext::element_markdown(size = 5),
    legend.text = element_text(angle = 33, vjust = 1, hjust = 1),
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing.x = unit(1, "lines")
  )

save_figure(
  output_dir, "umaps_celloracle_markov",
  plot = p1, width = 3.9, height = 1.40
)


p2 <- df2c %>%
  filter(factor %in% fois, is.na(branch)) %>%
  mutate(factor = sprintf("*%s* KO", factor)) %>%
  bind_rows(
    tibble(factor = "")
  ) %>%
  ggplot(aes(x, y)) +
  geom_point(aes(color = score), stroke = 0, size = 0.6) +
  facet_wrap2(
    ~factor,
    axes = "all"
  ) +
  scale_color_continuous_divergingx(
    palette = "PiYG",
    breaks = c(-1, 0, 1),
    guide = guide_colorbar(
      theme = theme(
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(0.3, "lines"),
        legend.text = element_text(size = 5)
      )
    )
  ) +
  theme_dimred2(arrow = arrow) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    color = "PS"
  ) +
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = 5, vjust = 0.95),
    legend.position = "inside",
    legend.title.position = "left",
    legend.justification = c("left", "bottom"),
    legend.direction = "horizontal",
    legend.text.position = "bottom",
    legend.margin = margin(),
    strip.text = ggtext::element_markdown(size = 5),
    plot.margin = margin(0, 0, 0, 0)
  )


save_figure(
  output_dir, "umaps_celloracle_ps",
  plot = p2, width = 1.8, height = 1
)
