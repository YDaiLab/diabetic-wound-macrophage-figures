library(tidyverse)
source("assets/theme_nature.R")

theme_set(theme_nature())

# Revision R2-2 — alternative-method cross-check (panel h). PAGA (partition-based
# graph abstraction; src/revision/trajectory/run_paga.py) scores the connectivity
# of every cluster pair from the single-cell kNN graph, independent of the MST /
# pseudotime. If the joint trajectory is real, its MST edges should be the
# high-connectivity pairs under PAGA, in the joint data and within each condition.
# Shows PAGA connectivity of joint-MST edges vs all other cluster pairs, per group.
# ----------------------------------------------------------------------
output_dir <- "figures/revision/trajectory/supp_panels"
source("make_figures/cfg.R")

d <- read_csv("results/revision/paga_connectivity.csv", show_col_types = FALSE) %>%
  mutate(
    group = recode(group, joint = "Joint", wt = "Non-diabetic", db = "Diabetic") %>%
      factor(c("Joint", "Non-diabetic", "Diabetic")),
    kind = if_else(in_joint_mst, "joint MST edge", "other pair") %>%
      factor(c("joint MST edge", "other pair"))
  )
# ----------------------------------------------------------------------

mst_col <- "#0A574C"; oth_col <- "#B7C2BD"

set.seed(1)
p <- ggplot(d, aes(kind, connectivity)) +
  geom_boxplot(aes(color = kind), width = 0.55, outlier.shape = NA,
               linewidth = 0.35, fill = NA) +
  geom_jitter(aes(color = kind, size = kind), width = 0.16, stroke = 0, alpha = 0.7) +
  facet_wrap(~group) +
  scale_color_manual(values = c("joint MST edge" = mst_col, "other pair" = oth_col), guide = "none") +
  scale_size_manual(values = c("joint MST edge" = 0.9, "other pair" = 0.5), guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "PAGA connectivity") +
  theme(
    axis.title = element_text(size = fs_min),
    axis.text.y = element_text(size = fs_min),
    axis.text.x = element_text(size = fs_min, angle = 20, hjust = 1),
    strip.text = element_text(size = fs_small),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# 90 mm wide (3.543 in), 3 facets
save_figure(output_dir, "panel_H_paga", plot = p, width = 3.543, height = 1.9)
