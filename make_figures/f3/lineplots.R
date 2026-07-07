library(tidyverse)
source("assets/theme_nature.R")
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(ggnewscale)
library(ggarrow)
library(ggarchery)
library(ggrepel)

theme_set(theme_nature())

# ----------------------------------------------------------------------
output_dir <- "figures/f3"
source("make_figures/cfg.R")

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
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df()

infer_tree <-
  "box/results/lamian/infer_tree.rds" %>%
  readRDS()

df

dat <- branches %>%
  lapply(function(name) {
    as_tibble(infer_tree[["order"]][[name]]) %>%
      mutate(name = name)
  }) %>%
  bind_rows() %>%
  left_join(
    df %>%
      select(Row.names, pseudotime, condition, clusterid, sample),
    by = join_by(value == Row.names)
  )

dat_l <- dat %>%
  group_by(name, clusterid) %>%
  summarise_if(is.numeric, median)

p <- dat %>%
  ggplot(aes(x = pseudotime)) +
  geom_vline(
    aes(xintercept = pseudotime, color = clusterid),
    dat_l,
    show.legend = FALSE,
    lwd = 0.4,
    lty = 2
  ) +
  geom_text_repel(
    aes(label = clusterid, y = 1),
    dat_l,
    max.iter = 0,
    size = pt2mm(fs_small),
    color = "#000",
    bg.color = "#fff"
  ) +
  scale_color_manual(values = cluster_pal) +
  new_scale_color() +
  geom_density(
    aes(
      y = after_stat(count),
      color = sample
    ),
    lwd = 0.8,
    key_glyph = "rect"
  ) +
  facet_nested_wrap(
    ~ condition + branches_names[name],
    nrow = 2,
    axes = "all",
    scales = "free",
    strip = strip_nested(
      text_x = list(
        element_text(size = fs_base, hjust = 0.5),
        element_text(size = fs_base, hjust = 0.5, margin = margin(t = 3, b = 6))
      ),
      by_layer_x = TRUE
    )
  ) +
  scale_color_manual(
    values = sample_pal,
    guide = guide_legend(byrow = TRUE)
  ) +
  labs(
    x = "Pseudotime",
    y = "Count-scaled density",
    color = ""
  ) +
  theme(
    axis.text = element_text(size = fs_small),
    legend.position = "top",
    legend.justification = "left",
    legend.location = "plot",
    legend.text.position = "right",
    legend.text = element_text(size = fs_base),
    legend.margin = margin(),
    legend.key.spacing = rel(1),
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2.2, "mm"),
    panel.spacing.y = unit(0.6, "lines")
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3), labels = scales::comma_format()) +
  coord_cartesian(clip = "off") +
  theme(
    ggh4x.facet.nestline = element_line(linewidth = 0.1)
  )

save_figure(output_dir, "lineplots", plot = p, width = 4, height = 3.2)
