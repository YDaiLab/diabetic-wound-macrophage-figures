library(tidyverse)
source("assets/theme_nature.R")
library(igraph)
library(colorspace)

theme_set(theme_nature())

# Revision R2-2 — schematic for panel A: how the best-match cell-Jaccard is
# computed. A real [joint edge x draw edge] Jaccard matrix for one representative
# draw (ND draw 3). Rows = joint MST edges (joint cluster labels), columns = that
# draw's MST edges (the draw's OWN, different cluster labels). Cell = cell-Jaccard
# of the two edges' cell sets. The row-max (best match) is outlined — that maximum
# is the "best-match cell-Jaccard" reported for each joint edge in panel F.
# ----------------------------------------------------------------------
output_dir <- "figures/revision/trajectory/supp_panels"
source("make_figures/cfg.R")

COND <- "wt"
DRAW <- 3 # hero draw
joint_file <- "box/results/lamian/infer_tree.rds"
draw_file <- sprintf("box/results/revision/lamian_condpca_draws/%s/seed%d/infer_tree.rds", COND, DRAW)
# ----------------------------------------------------------------------

jaccard <- function(A, B) length(intersect(A, B)) / length(union(A, B))
edge_sets <- function(res) {
  cl <- setNames(as.character(res$clusterid), names(res$clusterid))
  el <- igraph::as_edgelist(res$MSTtree)
  # order rows/cols by pseudotime-ish (min cluster id) so the diagonal-ish
  # best-match structure is legible
  keys <- apply(el, 1, function(r) paste(sort(as.integer(r)), collapse = "-"))
  sets <- lapply(seq_len(nrow(el)), function(i) names(cl)[cl %in% el[i, ]])
  names(sets) <- keys
  sets[order(keys)]
}

joint <- readRDS(joint_file)
draw <- readRDS(draw_file)

# restrict joint edges to this condition's cells (as the metric does)
cond_cells <- names(draw$clusterid)
jsets <- edge_sets(joint) %>% map(~ intersect(.x, cond_cells))
dsets <- edge_sets(draw)

mat <- expand_grid(joint = names(jsets), draw = names(dsets)) %>%
  mutate(jacc = map2_dbl(joint, draw, ~ jaccard(jsets[[.x]], dsets[[.y]]))) %>%
  group_by(joint) %>%
  mutate(is_best = jacc == max(jacc)) %>%
  ungroup() %>%
  mutate(
    joint = fct_rev(fct_inorder(joint)), # first joint edge at top
    draw = fct_inorder(draw)
  )

p <- ggplot(mat, aes(draw, joint)) +
  geom_tile(aes(fill = jacc), color = "white", linewidth = 0.2) +
  geom_tile(data = filter(mat, is_best), color = "#0A574C", fill = NA, linewidth = 0.5) +
  colorspace::scale_fill_continuous_sequential(
    palette = "Grays", rev = TRUE, begin = 0.05, end = 0.95, limits = c(0, 1),
    breaks = c(0, 1), labels = c("0", "1"), name = "cell-Jaccard",
    guide = guide_colorbar(theme = theme(
      legend.key.height = unit(0.3, "lines"), legend.key.width = unit(1.4, "lines"),
      legend.ticks = element_blank()
    ))
  ) +
  scale_x_discrete(position = "top") +
  coord_fixed(expand = FALSE) +
  labs(x = "run edges", y = "joint edges") +
  theme(
    axis.text.x = element_text(size = fs_min, angle = 90, hjust = 0, vjust = 0.5),
    axis.text.y = element_text(size = fs_min),
    axis.title = element_text(size = fs_min),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.location = "plot", # legend spans the whole plot width, not just panel
    legend.title = element_text(size = fs_min, vjust = 1),
    legend.text = element_text(size = fs_min),
    legend.margin = margin(t = -2)
  )

# 32 mm wide (1.26 in)
save_figure(output_dir, "schematic_jaccard_matrix", plot = p, width = 1.26, height = 1.43)
