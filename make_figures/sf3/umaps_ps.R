library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(ggrepel)
library(colorspace)
library(ggarrow)
library(ggarchery)
library(ggrastr)

theme_set(theme_nature())
options(ggrepel.max.overlaps = 6)
# ----------------------------------------------------------------------
output_dir <- "figures/sf3"
source("make_figures/cfg.R")
source("make_figures/load_edge_tbl.R")
source("make_figures/ggh4x_patch.R")

branches <- c(
  "branch: 12,6,9,3,4,8",
  "branch: 12,6,9,3,10"
)
rep_factors <- c("Cebpa", "Fli1", "Bhlhe40", "Kdm1a")
rep_factors <- fct_inorder(rep_factors)
rep_factors_pal <- setNames(
  RColorBrewer::brewer.pal(length(rep_factors), "Dark2"),
  rep_factors
)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df(pgd = TRUE)
tree <-
  "box/results/lamian/evaluate_uncertainty.rds" |>
  readRDS()

edge_tbl <- load_edge_tbl(tree, df) |>
  arrange(detection.rate)

df_ps <-
  "box/results/pgd_celloracle/perturb_scores.parquet" |>
  arrow::read_parquet() |>
  filter(
    is.na(branch),
    is.na(subbranch),
    factor %in% rep_factors
  ) |>
  mutate(
    factor = factor(factor, rep_factors)
  )

df2_ps <-
  "box/results/pgd_celloracle_overexpress/perturb_scores.parquet" |>
  arrow::read_parquet() |>
  filter(
    is.na(branch),
    is.na(subbranch),
    factor %in% rep_factors
  ) |>
  mutate(
    factor = factor(factor, rep_factors)
  )

p1 <- df_ps |>
  ggplot(aes(x, y)) +
  geom_point(aes(color = score), stroke = 0, size = 0.6) +
  scale_color_continuous_divergingx(
    palette = "PiYG",
    breaks = scales::pretty_breaks(n = 3),
    guide = guide_colorbar(
      theme = theme(
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(0.4, "lines")
      )
    ),
    limits = c(-1.05, 1.05),
    oob = scales::squish
  ) +
  facet_wrap2(
    ~factor,
    ncol = 2,
    axes = "all",
    strip = strip_themed(
      text_x = lapply(rep_factors, function(x) {
        ggtext::element_markdown(color = rep_factors_pal[x])
      })
    ),
    labeller = labeller(factor = function(x) sprintf("*%s* KO", x))
  ) +
  theme_embedding(arrow = arrow) +
  geom_arrow_segment(
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    edge_tbl,
    lwd = 0.3,
    position = position_attractsegment(
      start_shave = 0.2,
      end_shave = 0.2
    )
  ) +
  geom_shadowtext(
    aes(x = umap_1, y = umap_2, label = clusterid),
    df |>
      group_by(clusterid) |>
      summarise_if(is.numeric, median),
    color = "#000",
    bg.color = "#fff",
    size = pt2mm(fs_min)
  ) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    color = "Perturbation score (PS)"
  ) +
  theme(
    axis.title = element_text(size = fs_min),
    legend.title = element_text(size = fs_small, vjust = 0.95),
    legend.text = element_text(size = fs_small),
    legend.position = "bottom",
    legend.title.position = "left",
    legend.text.position = "bottom",
    legend.margin = margin(),
    strip.text = element_text(size = fs_small)
  )

p2 <- df2_ps |>
  ggplot(aes(x, y)) +
  geom_point(aes(color = score), stroke = 0, size = 0.6) +
  scale_color_continuous_divergingx(
    palette = "PiYG",
    breaks = scales::pretty_breaks(n = 3),
    guide = guide_colorbar(
      theme = theme(
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(0.4, "lines")
      )
    ),
    limits = c(-1.05, 1.05),
    oob = scales::squish
  ) +
  facet_wrap2(
    ~factor,
    ncol = 2,
    axes = "all",
    strip = strip_themed(
      text_x = lapply(rep_factors, function(x) {
        ggtext::element_markdown(color = rep_factors_pal[x])
      })
    ),
    labeller = labeller(factor = function(x) sprintf("*%s* overexpression", x))
  ) +
  theme_embedding(arrow = arrow) +
  geom_arrow_segment(
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    edge_tbl,
    lwd = 0.3,
    position = position_attractsegment(
      start_shave = 0.2,
      end_shave = 0.2
    )
  ) +
  geom_shadowtext(
    aes(x = umap_1, y = umap_2, label = clusterid),
    df |>
      group_by(clusterid) |>
      summarise_if(is.numeric, median),
    color = "#000",
    bg.color = "#fff",
    size = pt2mm(fs_min)
  ) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    color = "Perturbation score (PS)"
  ) +
  theme(
    axis.title = element_text(size = fs_min),
    legend.title = element_text(size = fs_small, vjust = 0.95),
    legend.text = element_text(size = fs_small),
    legend.position = "bottom",
    legend.title.position = "left",
    legend.text.position = "bottom",
    legend.margin = margin(),
    strip.text = element_text(size = fs_small)
  )

save_figure(output_dir, "umaps_ps_ko", plot = p1, width = 2.5, height = 3.15)
save_figure(output_dir, "umaps_ps_overexpression", plot = p2, width = 2.5, height = 3.15)
