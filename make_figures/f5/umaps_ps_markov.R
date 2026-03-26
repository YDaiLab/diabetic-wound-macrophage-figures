library(tidyverse)
library(ggbrandon)
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(ggrepel)
library(colorspace)
library(ggarrow)
library(ggarchery)
library(ggrastr)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 7))
options(ggrepel.max.overlaps = 6)
# ----------------------------------------------------------------------
output_dir <- "figures/f5"
source("make_figures/cfg.R")
source("make_figures/load_edge_tbl.R")
source("make_figures/ggh4x_patch.R")

branches <- c(
  "branch: 12,6,7,5",
  "branch: 12,6,9,3,10"
)
rep_factors <- c("Cebpa", "Nr1h3", "Irf4", "Kdm1a")
rep_factors <- fct_inorder(rep_factors)
rep_factors_pal <- setNames(
  RColorBrewer::brewer.pal(length(rep_factors), "Dark2"),
  rep_factors
)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df(pgd = TRUE)
tree <-
  "box/results/lamian/evaluate_uncertainty.rds" %>%
  readRDS()

edge_tbl <- load_edge_tbl(tree, df) %>%
  arrange(detection.rate)

# --- Perturbation scores (p1) ---
df_ps <-
  "box/results/pgd_celloracle/perturb_scores.parquet" %>%
  arrow::read_parquet() %>%
  filter(
    is.na(branch),
    is.na(subbranch),
    factor %in% rep_factors
  ) %>%
  mutate(
    factor = factor(factor, rep_factors)
  )

p1 <- df_ps %>%
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
  theme_dimred2(arrow = arrow) +
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
    df %>%
      group_by(clusterid) %>%
      summarise_if(is.numeric, median),
    color = "#000",
    bg.color = "#fff",
    size = 5 / .pt
  ) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    color = "Perturbation score (PS)"
  ) +
  theme(
    axis.title = element_text(size = 5),
    legend.title = element_text(size = 6, vjust = 0.95),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    legend.title.position = "left",
    legend.text.position = "bottom",
    legend.margin = margin(),
    strip.text = element_text(size = 6)
  )

# --- Markov chain transitions (p2) ---
n_bins <- 18

df_0 <-
  "box/results/pgd_celloracle/mc_transitions/0.parquet" %>%
  arrow::read_parquet() %>%
  filter(factor %in% rep_factors) %>%
  mutate(factor = factor(factor, rep_factors))

df_500 <-
  "box/results/pgd_celloracle/mc_transitions/500.parquet" %>%
  arrow::read_parquet() %>%
  filter(factor %in% rep_factors) %>%
  mutate(factor = factor(factor, rep_factors))

# Shared bins across both timepoints
x_range <- range(c(df_0$umap_1, df_500$umap_1), na.rm = TRUE)
y_range <- range(c(df_0$umap_2, df_500$umap_2), na.rm = TRUE)
x_breaks <- seq(x_range[1], x_range[2], length.out = n_bins + 1)
y_breaks <- seq(y_range[1], y_range[2], length.out = n_bins + 1)

bin_counts <- function(d) {
  d %>%
    mutate(
      xbin = cut(umap_1, breaks = x_breaks, include.lowest = TRUE, labels = FALSE),
      ybin = cut(umap_2, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
    ) %>%
    count(factor, xbin, ybin, name = "n")
}

df_change <- bin_counts(df_500) %>%
  rename(n_after = n) %>%
  full_join(
    bin_counts(df_0) %>% rename(n_before = n),
    by = c("factor", "xbin", "ybin")
  ) %>%
  replace_na(list(n_after = 0, n_before = 0)) %>%
  mutate(
    delta = n_after - n_before,
    # Bin center coordinates
    x = (x_breaks[xbin] + x_breaks[xbin + 1]) / 2,
    y = (y_breaks[ybin] + y_breaks[ybin + 1]) / 2
  ) %>%
  filter(delta != 0)

asinh_scale <- 50
asinh_lim <- max(abs(asinh(df_change$delta / asinh_scale)))

p2 <- df_change %>%
  ggplot(aes(x, y, fill = asinh(delta / asinh_scale))) +
  geom_tile(
    width = diff(x_breaks)[1],
    height = diff(y_breaks)[1]
  ) +
  scale_fill_gradientn(
    colours = pals::coolwarm(100),
    limits = c(-asinh_lim, asinh_lim),
    breaks = scales::pretty_breaks(n = 3),
    guide = guide_colorbar(
      theme = theme(
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(0.4, "lines")
      )
    )
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
  theme_dimred2(arrow = arrow) +
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
    ),
    alpha = 1,
    inherit.aes = FALSE
  ) +
  geom_shadowtext(
    aes(x = umap_1, y = umap_2, label = clusterid),
    df %>%
      group_by(clusterid) %>%
      summarise_if(is.numeric, median),
    color = "#000",
    bg.color = "#fff",
    size = 5 / .pt,
    alpha = 1,
    inherit.aes = FALSE
  ) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    fill = expression(asinh(Delta ~ count / 50))
  ) +
  theme(
    axis.title = element_text(size = 5),
    legend.title = element_text(size = 6, vjust = 1.05),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    legend.title.position = "left",
    legend.text.position = "bottom",
    legend.margin = margin(),
    strip.text = element_text(size = 6)
  )

# --- Combine and save ---
combined <- p1 + plot_spacer() + p2 +
  plot_layout(nrow = 1, widths = c(1, 0.15, 1)) &
  theme(plot.margin = margin(0, 0, 0, 0))

save_figure(output_dir, "umaps_ps_markov", plot = combined, width = 5.05, height = 3.15)
