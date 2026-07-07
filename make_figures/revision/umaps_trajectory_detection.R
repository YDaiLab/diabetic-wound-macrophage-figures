library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(patchwork)
library(legendry)
library(shadowtext)
library(ggrepel)
library(ggnewscale)
library(ggarrow)
library(ggarchery)
library(ggrastr)
library(colorspace)

theme_set(theme_nature())

# Revision task T1 / R2-2 (editor-flagged) — Variant A figure: per-condition
# bootstrap DETECTION RATE for the joint-cluster MST (workflow/scripts/
# revision/infer_tree_fixedclusters.R), shown on the shared PGD layout.
#
# Unlike umaps_trajectory_conditions.R (variant C, which re-infers clusters
# per condition in a condition-specific PC space), variant A holds the joint
# clusters and joint 10-PC space FIXED, so cluster anchor points are IDENTICAL
# across ALL THREE columns — there is one edge_tbl per condition, not a
# re-embedding.
#
#   columns: Reference (joint ND+DB)  |  Non-diabetic (wt)  |  Diabetic (db)
#   single row: clusters + detection-shaded backbone
#
# All three columns are shaded/annotated by variant A's per-edge, fixed-cluster
# bootstrap detection rate (results/revision/lamian_fixedclusters/
# {joint,wt,db}/edge_detection.csv) — one procedure throughout, so the joint
# reference is directly comparable to the conditions. The condition columns
# include condition-specific rewirings (e.g. wt: 1-3, 5-11 replacing 1-4/1-11)
# alongside joint edges that fail in that condition.
# ----------------------------------------------------------------------
output_dir <- "figures/revision"
source("make_figures/cfg.R")

# All three columns use variant A's per-edge fixed-cluster bootstrap detection
# rate; the joint file is the same procedure run on all cells (not the older
# per-branch Lamian evaluate_uncertainty), so a/d are directly comparable to
# b/c/e/f.
edge_files <- c(
  joint = "box/results/revision/lamian_fixedclusters/joint/edge_detection.csv",
  wt = "box/results/revision/lamian_fixedclusters/wt/edge_detection.csv",
  db = "box/results/revision/lamian_fixedclusters/db/edge_detection.csv"
)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df(pgd = TRUE) # umap_1/umap_2, clusterid (joint), condition, pseudotime

# Cluster anchor points and pseudotime from ALL cells. Joint clusters are fixed
# across conditions by design, so every panel places nodes at the same
# coordinates and shares one canonical early->late arrow orientation.
ctr <- df %>%
  group_by(clusterid) %>%
  summarise(
    umap_1 = median(umap_1), umap_2 = median(umap_2),
    pt = median(pseudotime), .groups = "drop"
  )
map_x <- setNames(ctr$umap_1, as.character(ctr$clusterid))
map_y <- setNames(ctr$umap_2, as.character(ctr$clusterid))
cluster_pt <- setNames(ctr$pt, as.character(ctr$clusterid))

# Orient each edge early->late by cluster-median pseudotime and attach layout
# geometry. The MST edge keys in edge_detection.csv are numerically sorted
# ("a|b", a<b), so raw source/target carry NO trajectory direction; without
# this, half the condition arrows point backwards along pseudotime.
finalize_edges <- function(edges) {
  edges %>%
    mutate(source = as.character(source), target = as.character(target)) %>%
    mutate(
      rev = cluster_pt[source] > cluster_pt[target],
      source0 = source, target0 = target,
      source = if_else(rev, target0, source0),
      target = if_else(rev, source0, target0),
      x = map_x[source], y = map_y[source],
      xend = map_x[target], yend = map_y[target],
      # perpendicular offset so the % label sits beside the arrow, not on it
      dx = xend - x, dy = yend - y, len = sqrt(dx^2 + dy^2),
      mid_x = (x + xend) / 2 - (dy / len) * 0.9,
      mid_y = (y + yend) / 2 + (dx / len) * 0.9
    ) %>%
    select(-rev, -source0, -target0) %>%
    filter(!is.na(x), !is.na(xend))
}

edge_tbl_from_csv <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(pct_label = scales::percent(detection_rate, accuracy = 1)) %>%
    finalize_edges()
}

edge_tbls <- lapply(edge_files, edge_tbl_from_csv)

# One compact colorbar geometry, shared by the detection and pseudotime guides so
# both legends render at the same size (matches the pseudotime bar in panel d).
key_h <- unit(1.2, "lines")
key_w <- unit(0.35, "lines")
colorbar_theme <- theme(legend.key.height = key_h, legend.key.width = key_w)

# theme_embedding() returns list(theme, guides); combine with c() (NOT `+`,
# which collapses a list-plus-theme to NULL and silently drops the styling).
theme_panel <- function(show_legend) {
  c(
    theme_embedding(arrow = arrow),
    list(theme(
      plot.title = element_text(size = fs_base, hjust = 0.5),
      axis.title = element_text(size = fs_min),
      legend.position = if (show_legend) "inside" else "none",
      legend.justification = c(0.02, 1.02),
      legend.title = element_text(size = fs_small, hjust = 0.5, margin = margin(b = 5)),
      legend.text = element_text(size = fs_min),
      legend.margin = margin()
    ))
  )
}

detection_scale <- function(show_legend) {
  colorspace::scale_color_continuous_sequential(
    palette = "Grays", rev = TRUE, begin = 0.3,
    limits = c(0, 1),
    breaks = scales::pretty_breaks(n = 3),
    labels = scales::percent_format(accuracy = 1),
    name = "Detection rate",
    guide = if (show_legend) guide_colorbar(theme = colorbar_theme) else "none"
  )
}

labels_layer <- function() {
  geom_shadowtext(
    aes(umap_1, umap_2, label = clusterid), ctr,
    color = "#000", bg.color = "#fff", size = pt2mm(fs_min), inherit.aes = FALSE
  )
}

# Repel the % labels off each other: several edges share a node (e.g. 6), so a
# fixed offset collides (12->6 "99%" and 6->2 "100%" overlapping at cluster 6).
# Anchor at the offset midpoint; a short connector is hidden (min.segment.length).
pct_labels_layer <- function(edge_tbl) {
  geom_text_repel(
    aes(mid_x, mid_y, label = pct_label), edge_tbl,
    color = "grey20", bg.color = "#fff", bg.r = 0.15,
    size = pt2mm(fs_min - 1), inherit.aes = FALSE,
    seed = 42, box.padding = 0.15, point.padding = 0,
    min.segment.length = 0.4, segment.size = 0.2, segment.color = "grey55",
    max.overlaps = Inf
  )
}

backbone_detection <- function(edge_tbl) {
  geom_arrow_segment(
    aes(x = x, y = y, xend = xend, yend = yend, color = detection_rate),
    edge_tbl, lwd = 0.4,
    position = position_attractsegment(start_shave = 0.25, end_shave = 0.25)
  )
}

point_scale <- function(color_var, show_legend) {
  if (color_var == "clusterid") {
    list(scale_color_manual(values = cluster_pal, guide = "none"))
  } else {
    list(
      scale_color_distiller(
        palette = "YlGnBu", direction = 1,
        breaks = c(range(df$pseudotime), median(df$pseudotime)),
        labels = c("Early", "Late", "Middle"),
        guide = if (show_legend) guide_colorbar(theme = colorbar_theme) else "none"
      ),
      labs(color = "Pseudotime")
    )
  }
}

cond_label <- function(cond) c(wt = "Non-diabetic", db = "Diabetic")[[cond]]

# One panel for all six cells. The reference (cond = "joint") draws every cell;
# a condition greys out the other condition's cells. All panels carry the same
# fixed-cluster detection-rate arrows + % labels. Legends are shown once:
# detection rate only in panel a, pseudotime only in panel d.
traj_panel <- function(cond, title, color_var,
                       show_detection_legend = FALSE, show_pt_legend = FALSE) {
  edge_tbl <- edge_tbls[[cond]]
  fg <- if (cond == "joint") df else df %>% filter(condition == cond_label(cond))

  p <- ggplot(mapping = aes(umap_1, umap_2))
  if (cond != "joint") {
    bg <- df %>% filter(condition != cond_label(cond))
    p <- p + rasterise(geom_point(data = bg, color = "grey88", stroke = 0, size = 0.35), dpi = 900)
  }
  p +
    rasterise(geom_point(aes(color = .data[[color_var]]), fg, stroke = 0, size = 0.4), dpi = 900) +
    point_scale(color_var, show_pt_legend) +
    new_scale_color() +
    backbone_detection(edge_tbl) +
    detection_scale(show_detection_legend) +
    pct_labels_layer(edge_tbl) +
    labels_layer() +
    labs(title = title, x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
    theme_panel(show_detection_legend || show_pt_legend)
}

# Single row: joint / ND / DB clusters + detection-shaded backbone.
# Detection-rate legend only in panel a. (Pseudotime coloring was dropped — the
# detection arrows carry the signal; the joint pseudotime added nothing here.)
fig <- (
  traj_panel("joint", "Reference (joint ND+DB)", "clusterid", show_detection_legend = TRUE) +
    traj_panel("wt", "Non-diabetic", "clusterid") +
    traj_panel("db", "Diabetic", "clusterid")
) +
  plot_annotation(tag_levels = "a") # plot.tag styling inherited from theme_nature()

# 124 mm wide (4.882 in), height scaled to preserve panel proportions.
save_figure(output_dir, "umaps_trajectory_detection", plot = fig, width = 4.882, height = 1.848)
