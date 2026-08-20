library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(shadowtext)
library(ggarrow)
library(ggarchery)
library(ggrastr)

theme_set(theme_nature())

# Revision R2-2 — BIRDccNEST alternative-method cross-check (panel i). The
# abstracted trajectory: BIRDccNEST's oMST backbone (a maximum-weight spanning
# tree over the community flow graph, rooted at the youngest community by
# pseudotime) drawn as directed arrows over the same ForceAtlas2 spring layout as
# the cell graph (panel h). Communities below MIN_CELLS are excluded.
#
# Two colorings on the SAME layout/backbone (cf. joint-cluster panels b-f):
#   - by BIRDccNEST community: community centroids labelled as the backbone nodes
#   - by the original joint Lamian cluster (cluster_pal): backbone arrows drawn
#     over cluster-colored cells with Lamian centroids labelled, so the reader
#     sees the inferred trajectory relative to the familiar clusters. (Arrows
#     still connect COMMUNITY centroids, shown as small grey nodes.)
#
# Direction concordance (all retained oMST edges increase in pseudotime) is
# reported in the figure legend, not on the panel.
# ----------------------------------------------------------------------
output_dir <- "figures/revision/trajectory/supp_panels"
source("make_figures/cfg.R")

bdir <- "box/results/revision/birdccnest/joint_res1p5"
MIN_CELLS <- 100 # res1.5: keep 5 communities >=100 cells (drops 7=66, 1=50, 6=46)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

community_pal <- c(
  "0" = "#4E79A7", "1" = "#F28E2B", "2" = "#E15759", "3" = "#76B7B2",
  "4" = "#59A14F", "5" = "#EDC948", "6" = "#B07AA1", "7" = "#FF9DA7",
  "8" = "#9C755F", "9" = "#BAB0AC", "10" = "#86BCB6", "11" = "#D37295",
  "12" = "#A0CBE8"
)

layout <- read_csv("results/revision/birdccnest_spring_layout.csv", show_col_types = FALSE) %>%
  select(cell, umap_1, umap_2)
comm_raw <- read_csv(file.path(bdir, "cell_communities.csv"), show_col_types = FALSE) %>%
  transmute(cell, community = as.character(community))
big <- comm_raw %>% count(community) %>% filter(n >= MIN_CELLS) %>% pull(community)

# community coloring (all cells; minor communities greyed = NA)
comm_all <- comm_raw %>%
  mutate(grp = factor(if_else(community %in% big, community, NA_character_),
                      intersect(names(community_pal), big)))
# big-community cells only, for centroid/backbone node positions
comm_big <- comm_raw %>% filter(community %in% big) %>%
  mutate(grp = factor(community, intersect(names(community_pal), big)))
# original joint Lamian cluster coloring (same source/palette as panels b-f)
lamian <- read_csv("results/pgd_cells_meta.csv", show_col_types = FALSE) %>%
  transmute(cell = Row.names, grp = factor(as.character(clusterid), names(cluster_pal)))

# community centroids = backbone node positions + arrow anchors
ctr <- layout %>% inner_join(comm_big, by = "cell") %>%
  group_by(grp) %>% summarise(x = median(umap_1), y = median(umap_2), .groups = "drop")
mx <- setNames(ctr$x, as.character(ctr$grp))
my <- setNames(ctr$y, as.character(ctr$grp))
lam_ctr <- layout %>% inner_join(lamian, by = "cell") %>%
  filter(!is.na(grp)) %>%
  group_by(grp) %>% summarise(x = median(umap_1), y = median(umap_2), .groups = "drop")

# oMST backbone (rooted parent -> child); edges to dropped communities fall out
omst <- read_csv(file.path(bdir, "omst_rooted_edges.csv"), show_col_types = FALSE) %>%
  mutate(x = mx[as.character(u)], y = my[as.character(u)],
         xend = mx[as.character(v)], yend = my[as.character(v)]) %>%
  filter(!is.na(x), !is.na(xend))

# oMST over cells colored by `grp`. comm_labeled = draw labelled community nodes
# (community coloring); else small grey backbone nodes + separate cluster labels.
omst_panel <- function(point_grp, pal, label_df, comm_labeled, title, out) {
  d <- layout %>% inner_join(point_grp, by = "cell")
  node_layers <- if (comm_labeled) {
    list(
      geom_point(aes(x, y), ctr, color = "#fff", size = 3.2, inherit.aes = FALSE),
      geom_point(aes(x, y, color = grp), ctr, size = 2.4, inherit.aes = FALSE),
      geom_shadowtext(aes(x, y, label = grp), ctr, color = "#000", bg.color = "#fff",
                      size = pt2mm(fs_min), inherit.aes = FALSE)
    )
  } else {
    list(geom_point(aes(x, y), ctr, color = "grey15", size = 1.8, inherit.aes = FALSE))
  }
  lab_layer <- if (!is.null(label_df)) {
    list(geom_shadowtext(aes(x, y, label = grp), label_df, color = "#000",
                         bg.color = "#fff", size = pt2mm(fs_min), inherit.aes = FALSE))
  } else {
    list()
  }
  p <- ggplot(d, aes(umap_1, umap_2)) +
    rasterise(geom_point(aes(color = grp), stroke = 0, size = 0.4, alpha = 0.3), dpi = 900) +
    scale_color_manual(values = pal, na.value = "grey85", guide = "none") +
    geom_arrow_segment(aes(x, y, xend = xend, yend = yend), omst, lwd = 0.5,
                       color = "grey10", inherit.aes = FALSE,
                       position = position_attractsegment(start_shave = 0.2, end_shave = 0.2)) +
    node_layers + lab_layer +
    labs(title = title, x = expression("FA"[1]), y = expression("FA"[2])) +
    c(theme_embedding(arrow = arrow), list(theme(
      plot.title = element_text(size = fs_base, hjust = 0.5),
      axis.title = element_text(size = fs_min), legend.position = "none"
    )))
  # ~33 mm square (1.28 in; 35% smaller than 1.97). Direction concordance reported in the legend.
  save_figure(output_dir, out, plot = p, width = 1.28, height = 1.28)
}

omst_panel(comm_all, community_pal, NULL, TRUE, "BIRDccNEST trajectory",
           "panel_I_birdccnest_omst")
omst_panel(lamian, cluster_pal, lam_ctr, FALSE, "Lamian clusters",
           "panel_I_birdccnest_omst_lamian")
