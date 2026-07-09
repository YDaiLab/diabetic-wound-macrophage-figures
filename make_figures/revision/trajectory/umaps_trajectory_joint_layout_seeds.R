library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(patchwork)
library(shadowtext)
library(ggarrow)
library(ggarchery)
library(ggrastr)

theme_set(theme_nature())

# Revision R2-2 — draw sweep. The condition-specific (variant C) trajectory was
# re-inferred under 10 raw k-means draws per condition
# (box/results/revision/lamian_condpca_seeds_raw/{wt,db}/draw{1..10}). Draw the
# shared-joint-layout view once PER DRAW so a good-looking consensus draw can be
# picked by eye. Each panel is annotated with the draw's cluster count (k) and
# its pseudotime concordance with the joint trajectory (rho, from raw_manifest.csv).
#   figures/revision/trajectory/draws/joint_layout_draw{N}.png
# ----------------------------------------------------------------------
output_dir <- "figures/revision/trajectory/draws"
source("make_figures/cfg.R")

draws <- 1:10
draw_dir <- "box/results/revision/lamian_condpca_seeds_raw"
manifest <- read_csv(file.path(draw_dir, "raw_manifest.csv"), show_col_types = FALSE)
cond_meta <- tribble(
  ~cond, ~condition,      ~label,
  "wt",  "Non-diabetic",  "Non-diabetic",
  "db",  "Diabetic",      "Diabetic"
)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df(pgd = TRUE) # joint PGD layout + condition + Row.names

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

cond_panel <- function(condition, label, tree, subtitle) {
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
    rasterise(geom_point(data = bg, color = "grey88", stroke = 0, size = 0.35), dpi = 300) +
    rasterise(geom_point(aes(color = factor(node)), fg, stroke = 0, size = 0.4), dpi = 300) +
    scale_color_manual(values = pal, na.value = "grey85", guide = "none") +
    backbone(edge_tbl(res, cn, fg_layout)) +
    labels_layer(node_labels(cn, fg_layout)) +
    labs(title = label, subtitle = subtitle,
         x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
    theme_panel
}

for (d in draws) {
  panels <- pmap(cond_meta, function(cond, condition, label) {
    m <- manifest %>% filter(condition == cond, draw == d)
    subtitle <- sprintf("k = %d   ρ = %.2f", m$n_clusters[1], m$rho_joint[1])
    cond_panel(condition, label,
               file.path(draw_dir, cond, paste0("draw", d), "infer_tree.rds"),
               subtitle)
  })
  fig <- (panels[[1]] + panels[[2]]) +
    plot_annotation(
      title = sprintf("draw %d", d), tag_levels = "a",
      theme = theme(plot.title = element_text(size = fs_base, face = "bold", hjust = 0.5))
    )
  ggsave(file.path(output_dir, sprintf("joint_layout_draw%02d.png", d)),
         fig, width = 4.882, height = 2.8, dpi = 300)
  message("draw ", d, " done")
}
