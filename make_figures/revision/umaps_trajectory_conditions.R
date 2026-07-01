library(tidyverse)
library(ggbrandon)
library(patchwork)
library(shadowtext)
library(ggarrow)
library(ggarchery)
library(ggrastr)
library(colorspace)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 7))

# Revision task T1 / R2-2 (editor-flagged) — per-condition trajectory on the
# shared PGD layout, to show the joint ND+DB trajectory reproduces when each
# condition is inferred independently.
#
#   Reference (joint ND+DB)  |  Non-diabetic (wt)  |  Diabetic (db)
#   row 1: clusters + branch backbone
#   row 2: pseudotime + branch backbone
#
# All three columns share the fixed PGD layout (umap_1/umap_2 from load_df), so
# occupancy/pseudotime are visually comparable. The wt/db columns use each
# condition's OWN independently-inferred clusters, pseudotime and branches
# (self-consistent per condition) — the quantitative cross-condition agreement is
# in src/revision/trajectory_concordance.R. Cells of the other condition are
# drawn as light-grey background so the shared canvas stays legible.
#
# NOTE (design choice to confirm): condition panels colour by the per-condition
# k-means clusters. For direct colour comparability instead, colour the condition
# panels by the joint `clusterid` (cluster_pal) and drop the per-condition labels.
# ----------------------------------------------------------------------
output_dir <- "figures/revision"
source("make_figures/cfg.R")
source("make_figures/load_edge_tbl.R")

joint_tree_file <- "box/results/lamian/evaluate_uncertainty.rds"
cond_meta <- tribble(
  ~cond, ~label,             ~condition,       ~tree_file,                                           ~res_file,
  "wt",  "Non-diabetic",     "Non-diabetic",   "box/results/revision/lamian/wt/evaluate_uncertainty.rds", "box/results/revision/lamian/wt/infer_tree.rds",
  "db",  "Diabetic",         "Diabetic",       "box/results/revision/lamian/db/evaluate_uncertainty.rds", "box/results/revision/lamian/db/infer_tree.rds"
)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df(pgd = TRUE) # umap_1/umap_2, clusterid (joint), condition, pseudotime, Row.names

median_xy <- function(d, group) {
  d %>%
    group_by(clusterid = .data[[group]]) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
}

theme_panel <- theme_dimred2(arrow = arrow) +
  theme(
    plot.title = element_text(size = 7, hjust = 0.5),
    axis.title = element_text(size = 5),
    legend.position = "none"
  )

backbone <- function(edges) {
  geom_arrow_segment(
    aes(x = x, y = y, xend = xend, yend = yend), edges, lwd = 0.3,
    position = position_attractsegment(start_shave = 0.2, end_shave = 0.2)
  )
}

labels_layer <- function(label_df) {
  geom_shadowtext(
    aes(umap_1, umap_2, label = clusterid), label_df,
    color = "#000", bg.color = "#fff", size = 5 / .pt
  )
}

# ---- reference (joint) ------------------------------------------------------
joint_tree <- readRDS(joint_tree_file)
edge_joint <- load_edge_tbl(joint_tree, df)

p_ref_cl <- ggplot(df, aes(umap_1, umap_2)) +
  rasterise(geom_point(aes(color = clusterid), stroke = 0, size = 0.4), dpi = 900) +
  scale_color_manual(values = cluster_pal) +
  backbone(edge_joint) + labels_layer(median_xy(df, "clusterid")) +
  labs(title = "Reference (joint ND+DB)", x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
  theme_panel

p_ref_pt <- ggplot(df, aes(umap_1, umap_2)) +
  rasterise(geom_point(aes(color = pseudotime), stroke = 0, size = 0.4), dpi = 900) +
  scale_color_viridis_c(guide = "none") +
  backbone(edge_joint) +
  labs(title = "Reference (joint ND+DB)", x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
  theme_panel

# ---- per-condition ----------------------------------------------------------
cond_panels <- pmap(cond_meta, function(cond, label, condition, tree_file, res_file) {
  res <- readRDS(res_file)
  tree <- readRDS(tree_file)

  cl <- enframe(res$clusterid, name = "Row.names", value = "cl_cond") %>%
    mutate(cl_cond = fct_inseq(as.character(cl_cond)))
  pt <- enframe(res$pseudotime, name = "Row.names", value = "pt_cond")

  bg <- df %>% filter(condition != !!condition)
  fg <- df %>%
    filter(condition == !!condition) %>%
    left_join(cl, by = "Row.names") %>%
    left_join(pt, by = "Row.names")

  # per-condition backbone: reuse load_edge_tbl by making clusterid = cl_cond
  edges <- load_edge_tbl(tree, fg %>% mutate(clusterid = cl_cond))
  cond_pal <- setNames(
    scales::hue_pal()(nlevels(fg$cl_cond)), levels(fg$cl_cond)
  )

  bg_layer <- geom_point(data = bg, color = "grey88", stroke = 0, size = 0.35)

  p_cl <- ggplot(mapping = aes(umap_1, umap_2)) +
    rasterise(bg_layer, dpi = 900) +
    rasterise(geom_point(aes(color = cl_cond), fg, stroke = 0, size = 0.4), dpi = 900) +
    scale_color_manual(values = cond_pal) +
    backbone(edges) + labels_layer(median_xy(fg, "cl_cond")) +
    labs(title = label, x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
    theme_panel

  p_pt <- ggplot(mapping = aes(umap_1, umap_2)) +
    rasterise(bg_layer, dpi = 900) +
    rasterise(geom_point(aes(color = pt_cond), fg, stroke = 0, size = 0.4), dpi = 900) +
    scale_color_viridis_c(guide = "none") +
    backbone(edges) +
    labs(title = label, x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
    theme_panel

  list(cl = p_cl, pt = p_pt)
})

row_cl <- p_ref_cl + cond_panels[[1]]$cl + cond_panels[[2]]$cl
row_pt <- p_ref_pt + cond_panels[[1]]$pt + cond_panels[[2]]$pt

fig <- (row_cl / row_pt) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_figure(output_dir, "umaps_trajectory_conditions", plot = fig, width = 6.5, height = 4.6)
