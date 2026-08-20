library(tidyverse)
source("assets/theme_nature.R")

theme_set(theme_nature())

# Revision R2-5, panel b — sensitivity to the diffusion parameter alpha, which
# also carries the benchmark against a conventional UMAP: alpha = 0 is the
# undiffused case, so the two-layout comparison is the endpoints of this sweep.
#
# Computed on the same regenerated layouts shown in panel a, not on the published
# ones. That matters less here than it does for panel a: these metrics are
# structural rather than coordinate-specific and agree to <= 0.02 between the two
# lineages wherever they overlap (alpha 0 to 0.6), including the peak position.
#
# Two spaces per metric: the UMAP, which is what was published, and the diffused
# PCs, which are deterministic and isolate alpha from UMAP's stochasticity.
# Smoothness saturates by the published alpha while trustworthiness peaks near
# it, so 0.6 sits at the turn rather than at the edge of a monotone trend.
# ----------------------------------------------------------------------
input_file <- "results/revision/pgd_alpha_metrics.csv"
output_dir <- "figures/revision/pgd/supp_panels"
source("make_figures/cfg.R")

PUBLISHED_ALPHA <- 0.6
# Diffused PCs first: it is the deterministic quantity, and the UMAP curve is the
# same effect seen through the stochastic projection step
# dimensionality in the labels: it explains the offset between the curves, since
# trustworthiness scores the layout against the 50-dim PCA space it came from
space_pal <- c("Diffused PCs (50 dim)" = "#cab2d6", "UMAP (2 dim)" = "#6a3d9a")
# ----------------------------------------------------------------------

d <- read_csv(input_file, show_col_types = FALSE) %>%
  pivot_longer(
    c(pseudotime_smoothness, trustworthiness),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(
    space = factor(
      recode(space, "Diffused PCs" = "Diffused PCs (50 dim)", "UMAP" = "UMAP (2 dim)"),
      names(space_pal)
    ),
    metric = recode(metric,
      pseudotime_smoothness = "Pseudotime smoothness",
      trustworthiness = "Trustworthiness"
    )
  )

# each curve's alpha = 0 value, carried across as the no-PGD baseline: anything
# above its own dashed line is an improvement on not diffusing at all
baseline <- d %>%
  filter(alpha == 0) %>%
  select(metric, space, value)

fig <- ggplot(d, aes(alpha, value, color = space)) +
  geom_vline(xintercept = PUBLISHED_ALPHA, linewidth = 0.25, linetype = "22", color = "grey70") +
  geom_hline(
    aes(yintercept = value, color = space), baseline,
    linewidth = 0.25, linetype = "22", show.legend = FALSE
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  scale_color_manual(values = space_pal, name = NULL) +
  scale_x_continuous(breaks = c(0, 0.4, PUBLISHED_ALPHA, 0.9)) +
  facet_wrap(~metric, ncol = 1, scales = "free_y", axes = "all") +
  labs(x = expression(alpha), y = NULL) +
  theme(
    strip.text = element_text(size = fs_small),
    axis.title = element_text(size = fs_small),
    axis.text = element_text(size = fs_min),
    legend.text = element_text(size = fs_min),
    legend.key.height = unit(0.3, "cm"),
    legend.position = "bottom",
    legend.margin = margin(t = -4),
    panel.grid.minor = element_blank()
  )

# 60 mm x 80 mm: the narrow half of row 2, height matched to panel C
save_figure(output_dir, "panel_B_alpha_sweep", plot = fig, width = 2.3, height = 2.8)
