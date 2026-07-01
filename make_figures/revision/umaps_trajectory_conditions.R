library(tidyverse)
library(ggbrandon)
library(patchwork)
library(shadowtext)
library(ggarrow)
library(ggarchery)
library(ggrastr)
library(colorspace)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 7))

# Revision task T1 / R2-2 (editor-flagged) — Variant C figure: per-condition
# trajectory inferred in a CONDITION-SPECIFIC PC space, shown on the shared PGD
# layout so occupancy/pseudotime stay visually comparable.
#
#   Reference (joint ND+DB)  |  Non-diabetic (wt)  |  Diabetic (db)
#   row 1: clusters + branch backbone
#   row 2: pseudotime + branch backbone
#
# The wt/db columns use each condition's OWN independently-inferred clusters,
# pseudotime and branches; quantitative cross-condition agreement is in
# src/revision/trajectory_concordance.R. Cells of the other condition are drawn
# as light-grey background so the shared canvas stays legible.
#
# Edges are built from names(res$order) (branch cluster sequences), NOT from
# load_edge_tbl()'s parsing of detection.rate rownames — the latter silently
# yields ZERO edges for some trees (e.g. db, whose detection.rate has bare
# integer rownames instead of "c(a,b)" node strings).
# ----------------------------------------------------------------------
output_dir <- "figures/revision"
source("make_figures/cfg.R")

joint_res_file <- "box/results/lamian/infer_tree.rds"
cond_meta <- tribble(
  ~cond, ~label,          ~condition,      ~res_file,
  "wt",  "Non-diabetic",  "Non-diabetic",  "box/results/revision/lamian_condpca/wt/infer_tree.rds",
  "db",  "Diabetic",      "Diabetic",      "box/results/revision/lamian_condpca/db/infer_tree.rds"
)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df(pgd = TRUE) # umap_1/umap_2, clusterid (joint), condition, pseudotime, Row.names

median_xy <- function(d, group) {
  d %>%
    group_by(clusterid = .data[[group]]) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
}

# Robust backbone edges: parse the branch cluster sequences from names(res$order)
# and place each node at the median layout coord of its cells (`df$clusterid`
# must hold the same node labels the tree uses).
edge_tbl_from_order <- function(res, d) {
  seqs <- lapply(names(res$order), function(nm) as.integer(str_extract_all(nm, "[0-9]+")[[1]]))
  edges <- bind_rows(lapply(seqs, function(v) {
    if (length(v) >= 2) tibble(source = v[-length(v)], target = v[-1]) else NULL
  })) %>% distinct()
  ctr <- d %>%
    mutate(node = as.integer(as.character(clusterid))) %>%
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
joint_res <- readRDS(joint_res_file)
edge_joint <- edge_tbl_from_order(joint_res, df)

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

# ---- per-condition (condition-specific PCA) ---------------------------------
cond_panels <- pmap(cond_meta, function(cond, label, condition, res_file) {
  res <- readRDS(res_file)

  cl <- enframe(res$clusterid, name = "Row.names", value = "cl_cond") %>%
    mutate(cl_cond = fct_inseq(as.character(cl_cond)))
  pt <- enframe(res$pseudotime, name = "Row.names", value = "pt_cond") %>%
    distinct(Row.names, .keep_all = TRUE) # order() repeats shared backbone cells

  bg <- df %>% filter(condition != !!condition)
  fg <- df %>%
    filter(condition == !!condition) %>%
    left_join(cl, by = "Row.names") %>%
    left_join(pt, by = "Row.names")

  edges <- edge_tbl_from_order(res, fg %>% mutate(clusterid = cl_cond))
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
