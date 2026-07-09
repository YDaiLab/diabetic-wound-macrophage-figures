library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(shadowtext)
library(ggarrow)
library(ggarchery)
library(ggrastr)

theme_set(theme_nature())

# Revision R2-2 supplementary — panels B, C, D as separate files for Illustrator
# assembly. All colored/labeled by the JOINT clusters (one shared color key).
#   B: reference — joint PGD layout, joint clusters, joint trajectory
#   C: Non-diabetic draw 3 — its own PGD layout, joint colors, DRAW-3 trajectory
#   D: Diabetic  draw 3 — same, for DB
# ----------------------------------------------------------------------
output_dir <- "figures/revision/trajectory/supp_panels"
source("make_figures/cfg.R")

DRAW <- 3
draw_dir <- "box/results/revision/lamian_condpca_draws"
pgd_dir <- "results/revision/pgd_draws"

# joint cluster per cell + the joint PGD layout (for panel B and coloring)
meta <- read_csv("results/pgd_cells_meta.csv", show_col_types = FALSE)
joint_cl <- meta %>% transmute(cell = Row.names, node = as.integer(as.character(clusterid)))
joint_layout <- meta %>% transmute(cell = Row.names, umap_1, umap_2)
joint_res <- readRDS("box/results/lamian/infer_tree.rds")
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

edge_tbl <- function(res, cell_node, layout) {
  seqs <- lapply(names(res$order), function(nm) as.integer(str_extract_all(nm, "[0-9]+")[[1]]))
  edges <- bind_rows(lapply(seqs, function(v) {
    if (length(v) >= 2) tibble(source = v[-length(v)], target = v[-1]) else NULL
  })) %>% distinct()
  ctr <- cell_node %>%
    inner_join(layout, by = "cell") %>%
    group_by(node) %>%
    summarise(x = median(umap_1), y = median(umap_2), .groups = "drop")
  mx <- setNames(ctr$x, ctr$node); my <- setNames(ctr$y, ctr$node)
  edges %>%
    mutate(x = mx[as.character(source)], y = my[as.character(source)],
           xend = mx[as.character(target)], yend = my[as.character(target)]) %>%
    filter(!is.na(x), !is.na(xend))
}

node_labels <- function(cell_node, layout) {
  cell_node %>% inner_join(layout, by = "cell") %>%
    group_by(clusterid = node) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
}

base_panel <- function(layout, edges, label_df, cn_joint, title, subtitle = NULL, pt_size = 0.4) {
  d <- layout %>% left_join(cn_joint, by = "cell")
  ggplot(d, aes(umap_1, umap_2)) +
    rasterise(geom_point(aes(color = factor(node)), stroke = 0, size = pt_size), dpi = 900) +
    scale_color_manual(values = cluster_pal, na.value = "grey85", guide = "none") +
    geom_arrow_segment(aes(x = x, y = y, xend = xend, yend = yend), edges, lwd = 0.3,
                       position = position_attractsegment(start_shave = 0.2, end_shave = 0.2)) +
    geom_shadowtext(aes(umap_1, umap_2, label = clusterid), label_df,
                    color = "#000", bg.color = "#fff", size = pt2mm(fs_min)) +
    labs(title = title, x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
    c(theme_embedding(arrow = arrow), list(theme(
      plot.title = element_text(size = fs_base, hjust = 0.5),
      axis.title = element_text(size = fs_min), legend.position = "none"
    )))
}

# ---- Panel B: reference (joint) --------------------------------------------
pB <- base_panel(
  joint_layout, edge_tbl(joint_res, joint_cl, joint_layout),
  node_labels(joint_cl, joint_layout), joint_cl,
  "Joint trajectory"
)
save_figure(output_dir, "panel_B_reference", plot = pB, width = 1.575, height = 1.7)

# ---- Panels C, D: draw-3 per condition -------------------------------------
cond_meta <- tribble(
  ~cond, ~label,                     ~out,
  "wt",  "Non-diabetic trajectory",  "panel_C_nd_draw3",
  "db",  "Diabetic trajectory",      "panel_D_db_draw3"
)
pwalk(cond_meta, function(cond, label, out) {
  layout <- read_csv(sprintf("%s/%s_draw%02d_layout.csv", pgd_dir, cond, DRAW),
                     show_col_types = FALSE) %>% select(cell, umap_1, umap_2)
  res <- readRDS(file.path(draw_dir, cond, paste0("seed", DRAW), "infer_tree.rds"))
  cn_draw <- enframe(res$clusterid, name = "cell", value = "node") %>%
    mutate(node = as.integer(as.character(node)))
  cn_joint <- joint_cl %>% filter(cell %in% layout$cell)
  p <- base_panel(
    layout, edge_tbl(res, cn_draw, layout), node_labels(cn_joint, layout), cn_joint,
    label, pt_size = 0.8  # ~half as many cells as B -> 2x point size
  )
  save_figure(output_dir, out, plot = p, width = 1.575, height = 1.7)
})
