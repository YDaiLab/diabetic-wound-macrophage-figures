library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(patchwork)
library(shadowtext)
library(ggarrow)
library(ggarchery)
library(ggrastr)

theme_set(theme_nature())

# Revision R2-2 — draw sweep on per-draw PGD layouts. For every raw draw
# (10 per condition), PGD is run on that draw's own trajectory + own full PCA
# (src/revision/run_pgd_draws.py) to give it its own layout. Each panel shows
# THAT draw's trajectory arrows on THAT draw's layout, but cells are colored and
# labeled by the REFERENCE JOINT clusters — so the coloring is identical across
# all draws and only the trajectory/layout varies, making draws comparable.
#
#   Two panels (ND draws, DB draws), each a 3-column grid of the 10 draws.
# ----------------------------------------------------------------------
output_dir <- "figures/revision/trajectory/supp_panels"
source("make_figures/cfg.R")

draw_dir <- "box/results/revision/lamian_condpca_draws"
pgd_dir <- "results/revision/pgd_draws"
manifest <- read_csv("box/results/revision/condpca_draw_manifest.csv", show_col_types = FALSE)

# joint cluster label per cell (the reference coloring/labeling, constant across draws)
joint_cl <- read_csv("results/pgd_cells_meta.csv", show_col_types = FALSE) %>%
  transmute(cell = Row.names, node = as.integer(as.character(clusterid)))

cond_meta <- tribble(
  ~cond, ~condition,      ~abbr,
  "wt",  "Non-diabetic",  "ND",
  "db",  "Diabetic",      "DB"
)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# draw trajectory arrows: draw-cluster centroids in the draw's PGD layout
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

# joint-cluster labels at joint-cluster centroids in the draw's PGD layout
node_labels <- function(cell_node, layout) {
  cell_node %>%
    inner_join(layout, by = "cell") %>%
    group_by(clusterid = node) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
}

panel <- function(cond, condition, abbr, d) {
  layout <- read_csv(sprintf("%s/%s_draw%02d_layout.csv", pgd_dir, cond, d),
                     show_col_types = FALSE) %>% select(cell, umap_1, umap_2)
  res <- readRDS(file.path(draw_dir, cond, paste0("seed", d), "infer_tree.rds"))

  cn_draw <- enframe(res$clusterid, name = "cell", value = "node") %>%
    mutate(node = as.integer(as.character(node)))
  cn_joint <- joint_cl %>% filter(cell %in% layout$cell)
  d_joint <- layout %>% left_join(cn_joint, by = "cell")

  ggplot(d_joint, aes(umap_1, umap_2)) +
    rasterise(geom_point(aes(color = factor(node)), stroke = 0, size = 0.25), dpi = 200) +
    scale_color_manual(values = cluster_pal, na.value = "grey85", guide = "none") +
    backbone_draw(edge_tbl(res, cn_draw, layout)) +
    geom_shadowtext(aes(umap_1, umap_2, label = clusterid), node_labels(cn_joint, layout),
                    color = "#000", bg.color = "#fff", size = pt2mm(fs_min - 1)) +
    labs(title = sprintf("%s draw %d", abbr, d)) +
    c(theme_embedding(arrow = arrow), list(theme(
      plot.title = element_text(size = fs_min, hjust = 0.5),
      axis.title = element_blank(), axis.line = element_blank(),
      legend.position = "none"
    )))
}

backbone_draw <- function(edges) {
  geom_arrow_segment(
    aes(x = x, y = y, xend = xend, yend = yend), edges, lwd = 0.25,
    position = position_attractsegment(start_shave = 0.15, end_shave = 0.15)
  )
}

nd <- map(1:10, ~ panel("wt", "Non-diabetic", "ND", .x))
db <- map(1:10, ~ panel("db", "Diabetic", "DB", .x))

# Two separate panels (ND draws, DB draws), each 3 columns x 4 rows
# (10 draws -> last row has a single panel). 60 mm wide (2.362 in) each.
save_figure(output_dir, "panel_E_nd_draws", plot = wrap_plots(nd, ncol = 3),
            width = 2.362, height = 4.0)
save_figure(output_dir, "panel_E_db_draws", plot = wrap_plots(db, ncol = 3),
            width = 2.362, height = 4.0)
