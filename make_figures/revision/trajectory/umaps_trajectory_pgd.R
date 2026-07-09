library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(patchwork)
library(shadowtext)
library(ggarrow)
library(ggarchery)
library(ggrastr)

theme_set(theme_nature())

# Revision task T1 / R2-2 — Variant C, re-drawn to remove the "cells moved"
# confusion of the shared-layout version.
#
# Each condition gets its OWN PGD layout (pgdiffusion run independently on that
# condition's variant-C trajectory + its own PC space; see
# src/revision/run_pgd_condition.py). Then a 2x2:
#
#            Non-diabetic                    Diabetic
#   top  | a: own clusters + own traj | b: own clusters + own traj   |
#   bot  | c: JOINT clusters + JOINT  | d: JOINT clusters + JOINT    |
#          traj (same layout as a)       traj (same layout as b)
#
# Columns are conditions; rows are the SAME cells in the SAME condition-specific
# embedding, read two ways: top = the condition's own re-inference, bottom =
# re-colored by the original joint clusters with the joint trajectory overlaid.
# Because top and bottom share coordinates, the joint structure is visibly the
# same cells — nothing moved, only the labelling/trajectory changed.
# ----------------------------------------------------------------------
output_dir <- "figures/revision/trajectory"
source("make_figures/cfg.R")

joint_res <- readRDS("box/results/lamian/infer_tree.rds")
# joint cluster label per cell (Row.names = cell, clusterid = joint 1..12)
joint_cl <- read_csv("results/pgd_cells_meta.csv", show_col_types = FALSE) %>%
  transmute(cell = Row.names, node = as.integer(as.character(clusterid)))

cond_meta <- tribble(
  ~cond, ~label,          ~pgd,                                  ~tree,
  "wt",  "Non-diabetic",  "results/revision/pgd_wt_layout.csv",  "box/results/revision/lamian_condpca/wt/infer_tree.rds",
  "db",  "Diabetic",      "results/revision/pgd_db_layout.csv",  "box/results/revision/lamian_condpca/db/infer_tree.rds"
)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Backbone edges from names(res$order): the branch cluster sequences are already
# in pseudotime order (origin first), so source->target arrows flow early->late.
# Nodes placed at the median layout coord of each cluster's cells (cell_node =
# cell -> node; layout = cell -> umap_1/umap_2).
edge_tbl <- function(res, cell_node, layout) {
  seqs <- lapply(names(res$order), function(nm) as.integer(str_extract_all(nm, "[0-9]+")[[1]]))
  edges <- bind_rows(lapply(seqs, function(v) {
    if (length(v) >= 2) tibble(source = v[-length(v)], target = v[-1]) else NULL
  })) %>% distinct()
  ctr <- cell_node %>%
    inner_join(layout, by = "cell") %>%
    group_by(node) %>%
    summarise(x = median(umap_1), y = median(umap_2), .groups = "drop")
  mx <- setNames(ctr$x, ctr$node)
  my <- setNames(ctr$y, ctr$node)
  edges %>%
    mutate(
      x = mx[as.character(source)], y = my[as.character(source)],
      xend = mx[as.character(target)], yend = my[as.character(target)]
    ) %>%
    filter(!is.na(x), !is.na(xend))
}

node_labels <- function(cell_node, layout) {
  cell_node %>%
    inner_join(layout, by = "cell") %>%
    group_by(clusterid = node) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
}

colorbar_theme <- theme(legend.key.height = unit(1.2, "lines"),
                         legend.key.width = unit(0.35, "lines"))

theme_panel <- function(show_legend = FALSE) {
  c(
    theme_embedding(arrow = arrow),
    list(theme(
      plot.title = element_text(size = fs_base, hjust = 0.5),
      plot.subtitle = element_text(size = fs_min, hjust = 0.5, color = "grey40",
                                   margin = margin(b = 3)),
      axis.title = element_text(size = fs_min),
      legend.position = if (show_legend) "inside" else "none",
      legend.justification = c(0.98, 1.02),
      legend.title = element_text(size = fs_min, hjust = 0.5, margin = margin(b = 4)),
      legend.text = element_text(size = fs_min),
      legend.margin = margin()
    ))
  )
}

backbone <- function(edges) {
  geom_arrow_segment(
    aes(x = x, y = y, xend = xend, yend = yend), edges, lwd = 0.3,
    position = position_attractsegment(start_shave = 0.2, end_shave = 0.2)
  )
}

labels_layer <- function(label_df) {
  geom_shadowtext(
    aes(umap_1, umap_2, label = clusterid), label_df,
    color = "#000", bg.color = "#fff", size = pt2mm(fs_min)
  )
}

# bottom row: cells colored by (joint) cluster
cluster_panel <- function(d, edges, label_df, pal, title, subtitle) {
  ggplot(d, aes(umap_1, umap_2)) +
    rasterise(geom_point(aes(color = factor(node)), stroke = 0, size = 0.4), dpi = 900) +
    scale_color_manual(values = pal, na.value = "grey85", guide = "none") +
    backbone(edges) +
    labels_layer(label_df) +
    labs(title = title, subtitle = subtitle,
         x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
    theme_panel()
}

# top row: cells colored by condition-specific pseudotime (late-at-top legend
# shown once). Each panel maps its own pseudotime range Early -> Late.
pseudotime_panel <- function(d, edges, label_df, title, subtitle, show_legend) {
  ggplot(d, aes(umap_1, umap_2)) +
    rasterise(geom_point(aes(color = pt), stroke = 0, size = 0.4), dpi = 900) +
    scale_color_distiller(
      palette = "YlGnBu", direction = 1,
      breaks = c(range(d$pt, na.rm = TRUE), median(d$pt, na.rm = TRUE)),
      labels = c("Early", "Late", "Middle"), name = NULL,
      guide = if (show_legend) guide_colorbar(theme = colorbar_theme) else "none"
    ) +
    backbone(edges) +
    labels_layer(label_df) +
    labs(title = title, subtitle = subtitle,
         x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
    theme_panel(show_legend)
}

panels <- pmap(cond_meta, function(cond, label, pgd, tree) {
  layout <- read_csv(pgd, show_col_types = FALSE) %>% select(cell, umap_1, umap_2)
  res <- readRDS(tree)

  # cluster labels + node positions come from the condition's own tree; used to
  # place the (own) trajectory arrows even though cells are colored by pseudotime.
  cn_own <- enframe(res$clusterid, name = "cell", value = "node") %>%
    mutate(node = as.integer(as.character(node)))
  pt <- enframe(res$pseudotime, name = "cell", value = "pt") %>%
    distinct(cell, .keep_all = TRUE) # order() repeats shared backbone cells

  # --- top: condition-specific pseudotime + own trajectory ---
  d_own <- layout %>% left_join(pt, by = "cell")
  p_own <- pseudotime_panel(
    d_own, edge_tbl(res, cn_own, layout), node_labels(cn_own, layout),
    label, "condition-specific pseudotime + trajectory",
    show_legend = cond == "wt"
  )

  # --- bottom: joint clusters + joint trajectory, same layout ---
  cn_joint <- joint_cl %>% filter(cell %in% layout$cell)
  d_joint <- layout %>% left_join(cn_joint, by = "cell")
  p_joint <- cluster_panel(
    d_joint, edge_tbl(joint_res, cn_joint, layout), node_labels(cn_joint, layout),
    cluster_pal, label, "joint clusters + trajectory"
  )

  list(own = p_own, joint = p_joint)
})

fig <- (panels[[1]]$own + panels[[2]]$own) /
  (panels[[1]]$joint + panels[[2]]$joint) +
  plot_annotation(tag_levels = "a") # plot.tag styling inherited from theme_nature()

# 124 mm wide (4.882 in); 2 rows of square panels.
save_figure(output_dir, "umaps_trajectory_pgd", plot = fig, width = 4.882, height = 4.9)
