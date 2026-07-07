library(tidyverse)
source("assets/theme_nature.R")
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(colorspace)

theme_set(theme_nature())

# ----------------------------------------------------------------------
output_dir <- "figures/sf4"
source("make_figures/cfg.R")

TOP_N <- 5
SIG_CUTOFF <- 0.05

branches <- c(
  "branch: 12,6,7,5",
  "branch: 12,6,9,3,10",
  "branch: 12,6,9,3,4,8"
)
branches_names <- c(
  "branch 1: 12,6,7,5",
  "branch 2: 12,6,9,3,10",
  "branch 3: 12,6,9,3,4,8"
)
branches_names <- setNames(branches_names, branches)

time_groups <- c("Early", "Middle", "Late")
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv("results/lamian/rna_tde_go.csv") |>
  filter(
    branch_name %in% branches,
    `Adjusted P-value` < SIG_CUTOFF
  ) |>
  mutate(
    peak_time_group = fct(peak_time_group, time_groups),
    branch_label = fct(branches_names[branch_name], branches_names),
    neg_log10_p = -log10(`Adjusted P-value`),
    fill_val = neg_log10_p
  )

# Add dummy rows for any missing branch × peak_time_group combos
# so that empty facets still render
all_combos <- expand_grid(
  branch_name = branches,
  peak_time_group = time_groups
) |>
  mutate(
    branch_label = fct(branches_names[branch_name], branches_names),
    peak_time_group = fct(peak_time_group, time_groups)
  )

missing <- all_combos |>
  anti_join(df, by = c("branch_name", "peak_time_group")) |>
  mutate(
    Term = " ",
    `Odds Ratio` = 0,
    `Adjusted P-value` = NA_real_,
    fill_val = NA_real_,
    is_placeholder = TRUE
  )

df <- df |>
  mutate(is_placeholder = FALSE) |>
  bind_rows(missing)

p <- df |>
  ggplot(aes(`Odds Ratio`, Term)) +
  # Colored bars for significant terms
  geom_col(aes(fill = fill_val)) +
  # "No significant terms" label in empty panels
  geom_text(
    data = \(d) filter(d, is_placeholder),
    aes(label = "No significant terms", y = Term),
    x = 0.1,
    hjust = 0.5, size = pt2mm(fs_min), color = "grey50"
  ) +
  # Reference line at OR = 1 (no enrichment), only in non-placeholder panels
  geom_vline(
    aes(xintercept = 1),
    lty = 2,
    lwd = 0.25,
    data = \(d) filter(d, !is_placeholder) |>
      distinct(branch_label, peak_time_group)
  ) +
  facet_nested_wrap(
    ~ branch_label + peak_time_group,
    nrow = 3,
    axes = "all",
    scales = "free",
    strip = strip_nested(
      text_x = list(
        element_text(size = fs_base, hjust = 0.5),
        element_text(size = fs_base, hjust = 0.5, margin = margin(t = 5, b = 5))
      ),
      by_layer_x = TRUE
    )
  ) +
  theme(
    ggh4x.facet.nestline = element_line(linewidth = 0.1),
    legend.position = "bottom",
    axis.text.x = element_text(size = fs_min),
    axis.text.y = element_text(size = fs_small),
    plot.margin = margin(r = 0.5, unit = "lines"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = fs_small, vjust = 0.9),
  ) +
  scale_fill_continuous_sequential(
    na.value = "gray80",
    name = expression(-log[10] ~ "(q)"),
    palette = "Blues 3", rev = TRUE, begin = 0.3,
    breaks = scales::pretty_breaks(n = 3),
    guide = guide_colorbar(
      theme = theme(
        legend.key.width  = unit(2.5, "lines"),
        legend.key.height = unit(0.4, "lines")
      )
    )
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.2)),
    breaks = scales::pretty_breaks(n = 3)
  ) +
  labs(y = NULL)

save_figure(output_dir, "barplot_go_terms", plot = p, width = 7.1, height = 6.5)
