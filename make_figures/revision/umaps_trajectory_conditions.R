library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(patchwork)
library(shadowtext)
library(ggnewscale)
library(ggarrow)
library(ggarchery)
library(ggrastr)
library(colorspace)

theme_set(theme_nature())

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
# Backbone edges are colored by each tree's OWN bootstrap detection rate
# (Lamian::evaluate_uncertainty, re-kmeans resampling) — the per-condition trees
# and the joint tree each use their own condition-specific bootstrap, which is
# the natural robustness measure for this fully-independent re-inference. (The
# joint tree's 3-4 / 4-8 edges read ~50% because Lamian decomposes detection at
# a branch fork; this is a decomposition artifact of the re-kmeans bootstrap, not
# genuine instability.)
#
# Edges are built from names(res$order) (branch cluster sequences), NOT from
# load_edge_tbl()'s parsing of detection.rate rownames — the latter silently
# yields ZERO edges for some trees. Detection is then matched onto those edges
# by parsing the detection.rate rownames (same "c(a,b)" / "a:b" grammar).
# ----------------------------------------------------------------------
output_dir <- "figures/revision"
source("make_figures/cfg.R")

joint_res_file <- "box/results/lamian/infer_tree.rds"
joint_eval_file <- "box/results/lamian/evaluate_uncertainty.rds"
cond_meta <- tribble(
  ~cond, ~label,          ~condition,      ~res_file,                                              ~eval_file,
  "wt",  "Non-diabetic",  "Non-diabetic",  "box/results/revision/lamian_condpca/wt/infer_tree.rds", "box/results/revision/lamian_condpca/wt/evaluate_uncertainty.rds",
  "db",  "Diabetic",      "Diabetic",      "box/results/revision/lamian_condpca/db/infer_tree.rds", "box/results/revision/lamian_condpca/db/evaluate_uncertainty.rds"
)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df(pgd = TRUE) # umap_1/umap_2, clusterid (joint), condition, pseudotime, Row.names

median_xy <- function(d, group) {
  d %>%
    group_by(clusterid = .data[[group]]) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
}

# Per-edge detection from an evaluate_uncertainty result: parse each detection.rate
# rowname ("c(6, 9, 3)" or "3:4") into consecutive-node edges, keyed by the
# canonical (min|max) pair, so it can be matched onto directed backbone edges.
edge_detection <- function(eval_res) {
  dr <- eval_res$detection.rate
  vals <- setNames(dr$detection.rate, rownames(dr))
  parse_nodes <- function(x) {
    if (grepl("^c\\(", x)) as.integer(str_extract_all(x, "\\d+")[[1]]) else as.integer(str_split(x, ":")[[1]])
  }
  imap_dfr(vals, function(d, nm) {
    v <- parse_nodes(nm)
    if (length(v) < 2) return(NULL)
    tibble(a = v[-length(v)], b = v[-1], detection_rate = d)
  }) %>%
    mutate(key = paste(pmin(a, b), pmax(a, b), sep = "|")) %>%
    distinct(key, .keep_all = TRUE) %>%
    select(key, detection_rate)
}

# Robust backbone edges: parse the branch cluster sequences from names(res$order)
# and place each node at the median layout coord of its cells (`df$clusterid`
# must hold the same node labels the tree uses). Detection is joined by edge key.
edge_tbl_from_order <- function(res, d, det) {
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
      key = paste(pmin(source, target), pmax(source, target), sep = "|"),
      x = mx[as.character(source)], y = my[as.character(source)],
      xend = mx[as.character(target)], yend = my[as.character(target)]
    ) %>%
    left_join(det, by = "key") %>%
    filter(!is.na(x), !is.na(xend))
}

# Compact colorbar geometry (matches umaps_trajectory_detection.R). The key
# dimensions must live in the GUIDE's theme, not just the panel theme, or a bare
# guide_colorbar() renders oversized.
colorbar_theme <- theme(legend.key.height = unit(1.2, "lines"), legend.key.width = unit(0.35, "lines"))

theme_panel <- function(show_legend = FALSE) {
  c(
    theme_embedding(arrow = arrow),
    list(theme(
      plot.title = element_text(size = fs_base, hjust = 0.5),
      axis.title = element_text(size = fs_min),
      legend.position = if (show_legend) "inside" else "none",
      legend.justification = c(0.02, 1.02),
      legend.title = element_text(size = fs_min, hjust = 0.5, margin = margin(b = 4)),
      legend.text = element_text(size = fs_min),
      legend.margin = margin()
    ))
  )
}

backbone <- function(edges) {
  geom_arrow_segment(
    aes(x = x, y = y, xend = xend, yend = yend, color = detection_rate), edges, lwd = 0.3,
    position = position_attractsegment(start_shave = 0.2, end_shave = 0.2)
  )
}

detection_scale <- function(show_legend) {
  colorspace::scale_color_continuous_sequential(
    palette = "Grays", rev = TRUE, begin = 0.3, limits = c(0, 1),
    breaks = scales::pretty_breaks(n = 3), labels = scales::percent_format(accuracy = 1),
    name = "Detection rate", na.value = "grey70",
    guide = if (show_legend) guide_colorbar(theme = colorbar_theme) else "none"
  )
}

labels_layer <- function(label_df) {
  geom_shadowtext(
    aes(umap_1, umap_2, label = clusterid), label_df,
    color = "#000", bg.color = "#fff", size = pt2mm(fs_min)
  )
}

# ---- reference (joint) ------------------------------------------------------
joint_res <- readRDS(joint_res_file)
edge_joint <- edge_tbl_from_order(joint_res, df, edge_detection(readRDS(joint_eval_file)))

p_ref_cl <- ggplot(df, aes(umap_1, umap_2)) +
  rasterise(geom_point(aes(color = clusterid), stroke = 0, size = 0.4), dpi = 900) +
  scale_color_manual(values = cluster_pal, guide = "none") +
  new_scale_color() +
  backbone(edge_joint) + detection_scale(TRUE) + labels_layer(median_xy(df, "clusterid")) +
  labs(title = "Reference (joint ND+DB)", x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
  theme_panel(TRUE)

# ---- per-condition (condition-specific PCA) ---------------------------------
# Clusters colored by condition-specific k-means; backbone shaded by that tree's
# own bootstrap detection rate. (Pseudotime-colored panels were dropped — the
# clusters + detection arrows carry the comparison.)
cond_panels <- pmap(cond_meta, function(cond, label, condition, res_file, eval_file) {
  res <- readRDS(res_file)
  det <- edge_detection(readRDS(eval_file))

  cl <- enframe(res$clusterid, name = "Row.names", value = "cl_cond") %>%
    mutate(cl_cond = fct_inseq(as.character(cl_cond)))

  bg <- df %>% filter(condition != !!condition)
  fg <- df %>%
    filter(condition == !!condition) %>%
    left_join(cl, by = "Row.names")

  edges <- edge_tbl_from_order(res, fg %>% mutate(clusterid = cl_cond), det)
  cond_pal <- setNames(
    scales::hue_pal()(nlevels(fg$cl_cond)), levels(fg$cl_cond)
  )

  ggplot(mapping = aes(umap_1, umap_2)) +
    rasterise(geom_point(data = bg, color = "grey88", stroke = 0, size = 0.35), dpi = 900) +
    rasterise(geom_point(aes(color = cl_cond), fg, stroke = 0, size = 0.4), dpi = 900) +
    scale_color_manual(values = cond_pal, guide = "none") +
    new_scale_color() +
    backbone(edges) + detection_scale(FALSE) + labels_layer(median_xy(fg, "cl_cond")) +
    labs(title = label, x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
    theme_panel(FALSE)
})

fig <- (p_ref_cl + cond_panels[[1]] + cond_panels[[2]]) +
  plot_annotation(tag_levels = "a") # plot.tag styling inherited from theme_nature()

# 124 mm wide (4.882 in), single row of square panels.
save_figure(output_dir, "umaps_trajectory_conditions", plot = fig, width = 4.882, height = 1.848)
