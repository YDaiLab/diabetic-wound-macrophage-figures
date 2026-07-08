library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(patchwork)
library(shadowtext)
library(ggarrow)
library(ggarchery)
library(ggrastr)

theme_set(theme_nature())

# Revision task T1 / R2-2 — Variant C, "shared-layout" alternative to
# umaps_trajectory_pgd.R.
#
# Instead of giving each condition its own PGD embedding, this keeps the single
# JOINT PGD layout and just subsets it to ND / DB, colouring each condition's
# cells by its OWN condition-specific clusters and overlaying its OWN
# (condition-specific) trajectory. The other condition's cells are drawn light
# grey so the shared canvas stays legible.
#
#   Non-diabetic | Diabetic   (one row, joint PGD coordinates)
# ----------------------------------------------------------------------
output_dir <- "figures/revision"
source("make_figures/cfg.R")

df <- load_df(pgd = TRUE) # umap_1/umap_2 (JOINT PGD), clusterid (joint), condition, Row.names

cond_meta <- tribble(
  ~cond, ~condition,      ~label,          ~tree,
  "wt",  "Non-diabetic",  "Non-diabetic",  "box/results/revision/lamian_condpca/wt/infer_tree.rds",
  "db",  "Diabetic",      "Diabetic",      "box/results/revision/lamian_condpca/db/infer_tree.rds"
)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Condition-specific trajectory placed on the joint layout: branch cluster
# sequences from names(res$order) (already origin-first, so source->target flows
# early->late), nodes at the median joint-PGD coord of each cluster's cells.
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

theme_panel <- c(
  theme_embedding(arrow = arrow),
  list(theme(
    plot.title = element_text(size = fs_base, hjust = 0.5),
    plot.subtitle = element_text(size = fs_min, hjust = 0.5, color = "grey40",
                                 margin = margin(b = 3)),
    axis.title = element_text(size = fs_min),
    legend.position = "none"
  ))
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
    color = "#000", bg.color = "#fff", size = pt2mm(fs_min)
  )
}

panels <- pmap(cond_meta, function(cond, condition, label, tree) {
  res <- readRDS(tree)
  cn <- enframe(res$clusterid, name = "cell", value = "node") %>%
    mutate(node = as.integer(as.character(node)))

  bg <- df %>% filter(condition != !!condition)
  fg <- df %>%
    filter(condition == !!condition) %>%
    left_join(cn, by = c("Row.names" = "cell"))
  fg_layout <- fg %>% transmute(cell = Row.names, umap_1, umap_2)

  pal <- setNames(scales::hue_pal()(n_distinct(cn$node)), sort(unique(cn$node)))

  ggplot(mapping = aes(umap_1, umap_2)) +
    rasterise(geom_point(data = bg, color = "grey88", stroke = 0, size = 0.35), dpi = 900) +
    rasterise(geom_point(aes(color = factor(node)), fg, stroke = 0, size = 0.4), dpi = 900) +
    scale_color_manual(values = pal, na.value = "grey85", guide = "none") +
    backbone(edge_tbl(res, cn, fg_layout)) +
    labels_layer(node_labels(cn, fg_layout)) +
    labs(title = label, subtitle = "condition-specific clusters + trajectory",
         x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
    theme_panel
})

fig <- (panels[[1]] + panels[[2]]) +
  plot_annotation(tag_levels = "a") # plot.tag styling inherited from theme_nature()

# 124 mm wide (4.882 in); one row of square panels.
save_figure(output_dir, "umaps_trajectory_joint_layout", plot = fig, width = 4.882, height = 2.7)
