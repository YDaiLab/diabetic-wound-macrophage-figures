library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(ggrastr)

theme_set(theme_nature())

# Revision R2-5, panel a — the qualitative counterpart to the metrics in panel b:
# the layout at each alpha, coloured by pseudotime. alpha = 0 is the undiffused
# case and 0.6 is the value used in the paper.
#
# NOT the published layouts. The whole strip was regenerated in one environment
# (src/revision/pgd/alpha_sweep.py) so the six panels are mutually comparable.
# The diffused PCs are identical to the published ones, but the UMAP step is not:
# relative RMSD after Procrustes is 0.17 at alpha = 0.6 and 0.95 at alpha = 0,
# the latter because the published conventional UMAP came from Seurat's RunUMAP
# rather than scanpy. Same structure, different coordinates -- so label these
# "value used", never "published", and say so in the caption. The metrics in
# panel b are structural and agree to <= 0.02 between the two lineages.
#
# Coordinates are Procrustes-aligned onto the published layout. UMAP fixes
# neither orientation nor reflection, so without alignment each panel would
# arrive rotated at random and the progression would be unreadable. Alignment is
# display-only: it is a rigid transform, so it changes none of the metrics in
# panel b, and the layout used for the analysis is untouched.
# ----------------------------------------------------------------------
input_file <- "results/revision/pgd_alpha_sweep.csv"
output_dir <- "figures/revision/pgd/supp_panels"
source("make_figures/cfg.R")

PUBLISHED_ALPHA <- 0.6
# ----------------------------------------------------------------------

d <- read_csv(input_file, show_col_types = FALSE)

# strip labels carry the two reference points so the reader does not have to be
# told separately which panel is the published one
lab <- d %>%
  distinct(alpha) %>%
  arrange(alpha) %>%
  mutate(label = case_when(
    alpha == 0 ~ paste0("alpha == 0 ~ '(no diffusion)'"),
    alpha == PUBLISHED_ALPHA ~ paste0("alpha == ", alpha, " ~ '(value used)'"),
    TRUE ~ paste0("alpha == ", alpha)
  )) %>%
  deframe()

fig <- d %>%
  mutate(alpha = factor(lab[as.character(alpha)], lab)) %>%
  ggplot(aes(aligned_1, aligned_2, color = pseudotime)) +
  rasterise(geom_point(size = 0.22, stroke = 0), dpi = 900) +
  scale_color_distiller(
    palette = "YlGnBu", direction = 1,
    breaks = range(d$pseudotime), labels = c("Early", "Late"),
    guide = guide_colorbar(theme = theme(
      legend.key.height = unit(0.9, "lines"),
      legend.key.width = unit(0.3, "lines")
    ))
  ) +
  theme_embedding(arrow = arrow) +
  facet_wrap(~alpha, nrow = 1, labeller = label_parsed, axes = "all") +
  labs(color = "Pseudotime", x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
  theme(
    strip.text = element_text(size = fs_small),
    legend.title = element_text(size = fs_min),
    legend.text = element_text(size = fs_min),
    axis.title = element_text(size = fs_min),
    legend.margin = margin(l = -4),
    panel.spacing.x = unit(0.2, "lines")
  )

# 180 mm = 7.087 in full page width
save_figure(output_dir, "panel_A_interpolation", plot = fig, width = 7.087, height = 1.5)
