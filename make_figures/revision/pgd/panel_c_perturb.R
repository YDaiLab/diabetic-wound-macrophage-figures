library(tidyverse)
source("assets/theme_nature.R")

theme_set(theme_nature())

# Revision R2-5, panel c — do the biological conclusions depend on PGD?
#
# CellOracle was run twice on the same cells, the same GRN and the same
# perturbations, differing only in the embedding the vector field is projected
# onto. Each point is one transcriptional regulator, positioned by its rank for
# summed inhibitory perturbation score within a branch under each embedding;
# points on the diagonal are regulators the two runs agree about.
#
# RANKS, not magnitudes: CellOracle fits the vector-field grid to each embedding
# separately (8,598 vs 5,013 grid points per TR), so the sums are on different
# scales by construction. Spearman rho for both inhibitory and promoting scores
# is annotated per branch; the regulators named in the Results are filled.
# ----------------------------------------------------------------------
input_file <- "results/revision/pgd_perturb_comparison.csv"
output_dir <- "figures/revision/pgd/supp_panels"
source("make_figures/cfg.R")

# TRs the manuscript names as exemplars, per score type: Fig. 5a for the summed
# negative (inhibitory) PS, Fig. 5b for the summed positive (promoting) PS. The
# two sets differ, so filling the Fig. 5a set in both rows would mislabel the
# promotion row.
# NB the Results name "Flk1" among the strong promoters, but no such factor
# exists in the network or the CellOracle output; Fli1 does and ranks 4th
# (UMAP) / 1st (PGD) on branch 2. Read as Fli1 here.
# Same representative TRs and palette as make_figures/f5/scatterplots_ps.R, which
# highlights these four across both its negative- and positive-PS panels, so a
# reader can carry them straight from Fig. 5 into this panel.
rep_factors <- c("Cebpa", "Nr1h3", "Irf4", "Kdm1a")
rep_factors_pal <- setNames(
  RColorBrewer::brewer.pal(length(rep_factors), "Set2"), rep_factors
)
# Restricted to the three branches the manuscript interprets: Fig. 5a plots
# branch 1 vs 2 and Fig. S4 compares 2 vs 3. The remaining branch (12,6,2) and
# the backbone are not units of interpretation, so they are not shown here; the
# full five-branch comparison is in results/revision/pgd_perturb_comparison.csv.
branch_labels <- c(
  "branch: 12,6,7,5" = "Branch 1: 12-6-7-5",
  "branch: 12,6,9,3,10" = "Branch 2: 12-6-9-3-10",
  "branch: 12,6,9,3,4,8" = "Branch 3: 12-6-9-3-4-8"
)
# ----------------------------------------------------------------------

# Both score types, as the manuscript uses both: Fig. 5a is the summed negative
# (inhibitory) PS, Fig. 5b the summed positive (promoting) PS. Rank 1 = strongest
# in each case.
d <- read_csv(input_file, show_col_types = FALSE) %>%
  filter(perturbation == "Knockout", branch %in% names(branch_labels)) %>%
  mutate(branch = factor(branch_labels[branch], branch_labels)) %>%
  pivot_longer(
    starts_with("rank_"),
    names_to = c("kind", "layout"), names_pattern = "rank_(neg|pos)_(umap|pgd)",
    values_to = "rank"
  ) %>%
  pivot_wider(names_from = layout, values_from = rank) %>%
  mutate(
    kind = factor(
      # wording as Fig. 5's axis titles, which call the magnitude of the
      # summed negative scores the "negative PS sum"
      recode(kind,
        neg = "Inhibition (negative PS sum)",
        pos = "Promotion (positive PS sum)"
      ),
      c("Inhibition (negative PS sum)", "Promotion (positive PS sum)")
    ),
    rep = factor > "" & factor %in% rep_factors,
    # "" not NA: ggrepel only avoids points present in its own layer, and
    # na.rm drops them before repelling, so empty labels keep all 98 points
    # as obstacles while drawing no text (as f5/scatterplots_ps.R)
    label = if_else(rep, factor, "")
  ) %>%
  arrange(rep) # representative TRs drawn on top

# one rho per facet, in the top-left: points hug the diagonal, so it stays clear
rho <- d %>%
  group_by(branch, kind) %>%
  summarise(
    label = sprintf("rho=='%.2f'", cor(umap, pgd, method = "spearman")),
    .groups = "drop"
  )

fig <- ggplot(d, aes(umap, pgd)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.2, color = "grey80") +
  # point style as make_figures/f5/scatterplots_ps.R: grey fill, white outline,
  # translucent halo behind the representative TRs
  geom_point(
    data = ~ filter(.x, !rep),
    fill = "grey50", shape = 21, color = "white", size = 1.0, stroke = 0.2
  ) +
  geom_point(
    aes(color = rep_factors_pal[factor]),
    data = ~ filter(.x, rep),
    shape = 16, size = 3.5, alpha = 0.3
  ) +
  geom_point(
    aes(fill = rep_factors_pal[factor]),
    data = ~ filter(.x, rep),
    shape = 21, color = "white", size = 1.0, stroke = 0.2
  ) +
  ggrepel::geom_text_repel(
    aes(label = label),
    fontface = "italic",
    size = pt2mm(fs_min), segment.size = 0.25,
    box.padding = 0.3, point.padding = 0.2, # clearance around the 98 points
    min.segment.length = 0, max.overlaps = Inf, max.iter = 1e5, seed = 42
  ) +
  geom_text(
    aes(x = 3, y = Inf, label = label), rho,
    parse = TRUE, hjust = 0, vjust = 1.6, size = pt2mm(fs_min), inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(breaks = c(1, 50, 98)) +
  scale_y_continuous(breaks = c(1, 50, 98)) +
  coord_fixed(clip = "off") +
  facet_grid(kind ~ branch, switch = "y", axes = "all") +
  labs(
    x = "Perturbation score rank (1 = strongest), conventional UMAP",
    y = "Perturbation score rank (1 = strongest), PGD"
  ) +
  theme(
    strip.text = element_text(size = fs_small),
    strip.placement = "outside", # kind labels sit left of the axis, not against the panel
    strip.background = element_blank(),
    axis.title = element_text(size = fs_small),
    axis.text = element_text(size = fs_min),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.35, "lines") # keep adjacent facets' 98 / 1 ticks apart
  )

# 104 mm x 80 mm: the wide half of row 2 (coord_fixed, so height tracks width)
save_figure(output_dir, "panel_C_perturb", plot = fig, width = 4.5, height = 3.5)
