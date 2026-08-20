library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(shadowtext)
library(ggrastr)

theme_set(theme_nature())

# Revision R2-2 — BIRDccNEST alternative-method cross-check (panel h). BIRDccNEST
# is an independent community-flow trajectory method: it prunes the single-cell
# kNN graph, detects communities (Louvain), and extracts an oMST backbone. Run in
# the birdccnest workflow (box/results/revision/birdccnest/joint_res1p5).
#
# Panel h — the cell-cell graph itself: every cell is a node placed by a
# ForceAtlas2 spring layout of BIRDccNEST's own pruned graph
# (src/revision/trajectory/birdccnest_spring_layout.py), edges are that graph's
# strongest cell-cell connections. The abstracted backbone (panel i) and cluster
# correspondence (panel j) are separate panels.
#
# Two colorings on the SAME layout/edges (like the joint-cluster panels b-f):
#   - by BIRDccNEST community (minor communities <MIN_CELLS greyed)
#   - by the original joint Lamian cluster (cluster_pal, so it matches b-f)
# so the reader can anchor the independent-method graph to the familiar clusters.
# ----------------------------------------------------------------------
output_dir <- "figures/revision/trajectory/supp_panels"
source("make_figures/cfg.R")

bdir <- "box/results/revision/birdccnest/joint_res1p5"
MIN_CELLS <- 100 # res1.5: keep 5 communities >=100 cells (drops 7=66, 1=50, 6=46)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# qualitative community palette (Tableau20-ish), distinct from cluster_pal
community_pal <- c(
  "0" = "#4E79A7", "1" = "#F28E2B", "2" = "#E15759", "3" = "#76B7B2",
  "4" = "#59A14F", "5" = "#EDC948", "6" = "#B07AA1", "7" = "#FF9DA7",
  "8" = "#9C755F", "9" = "#BAB0AC", "10" = "#86BCB6", "11" = "#D37295",
  "12" = "#A0CBE8"
)

# spring layout + cell-cell edges (endpoints already in layout coords)
layout <- read_csv("results/revision/birdccnest_spring_layout.csv", show_col_types = FALSE) %>%
  select(cell, umap_1, umap_2)
edges <- read_csv("results/revision/birdccnest_spring_edges.csv", show_col_types = FALSE)

# all cells kept; big communities colored, minor ones greyed (community = NA)
comm_raw <- read_csv(file.path(bdir, "cell_communities.csv"), show_col_types = FALSE) %>%
  transmute(cell, community = as.character(community))
big <- comm_raw %>% count(community) %>% filter(n >= MIN_CELLS) %>% pull(community)
comm <- comm_raw %>%
  mutate(grp = factor(if_else(community %in% big, community, NA_character_),
                      intersect(names(community_pal), big)))

# original joint Lamian cluster per cell (same source/palette as panels b-f)
lamian <- read_csv("results/pgd_cells_meta.csv", show_col_types = FALSE) %>%
  transmute(cell = Row.names, grp = factor(as.character(clusterid), names(cluster_pal)))

# centroids (median layout position) per group, for the shadowtext labels
centroids <- function(d) {
  d %>% filter(!is.na(grp)) %>% group_by(grp) %>%
    summarise(x = median(umap_1), y = median(umap_2), .groups = "drop")
}

# cell graph colored by `grp`: shared edges + points + centroid labels
network_panel <- function(cell_grp, pal, title, out) {
  d <- layout %>% inner_join(cell_grp, by = "cell")
  p <- ggplot(d, aes(umap_1, umap_2)) +
    rasterise(geom_segment(aes(x, y, xend = xend, yend = yend), edges,
                           color = "grey65", linewidth = 0.05, alpha = 0.04,
                           inherit.aes = FALSE), dpi = 600) +
    rasterise(geom_point(aes(color = grp), stroke = 0, size = 0.4), dpi = 900) +
    scale_color_manual(values = pal, na.value = "grey85", guide = "none") +
    geom_shadowtext(aes(x, y, label = grp), centroids(d), color = "#000",
                    bg.color = "#fff", size = pt2mm(fs_min), inherit.aes = FALSE) +
    labs(title = title, x = expression("FA"[1]), y = expression("FA"[2])) +
    c(theme_embedding(arrow = arrow), list(theme(
      plot.title = element_text(size = fs_base, hjust = 0.5),
      axis.title = element_text(size = fs_min), legend.position = "none"
    )))
  # ~33 mm square (1.28 in; 35% smaller than 1.97)
  save_figure(output_dir, out, plot = p, width = 1.28, height = 1.28)
}

network_panel(comm, community_pal, "BIRDccNEST cell graph",
              "panel_H_birdccnest_network")
network_panel(lamian, cluster_pal, "Cell graph (Lamian clusters)",
              "panel_H_birdccnest_network_lamian")
