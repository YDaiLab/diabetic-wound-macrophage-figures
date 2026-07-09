library(tidyverse)
source("assets/theme_nature.R")
source("assets/theme_embedding.R")
library(patchwork)
library(shadowtext)
library(ggrastr)

theme_set(theme_nature())

# Revision R2-2 — PAGA graph (alternative-method cross-check). PAGA scores the
# connectivity of every cluster pair from the single-cell kNN graph, with no MST
# and no pseudotime (src/revision/trajectory/run_paga.py). Drawn as the classic
# node-link abstraction on the joint PGD layout: nodes = joint clusters at their
# layout centroids, edges weighted (width + opacity) by PAGA connectivity. Same
# node positions across Joint / ND / DB so the graphs are directly comparable;
# the joint MST backbone (panel b) lights up as the strong PAGA edges in all
# three. Faint edges are the low-connectivity pairs.
# ----------------------------------------------------------------------
output_dir <- "figures/revision/trajectory/supp_panels"
source("make_figures/cfg.R")

df <- load_df(pgd = TRUE) # joint PGD layout + joint clusterid
conn <- read_csv("results/revision/paga_connectivity.csv", show_col_types = FALSE)
# ----------------------------------------------------------------------

# node positions = joint-cluster centroids on the joint layout (shared by all groups)
ctr <- df %>%
  group_by(node = as.integer(as.character(clusterid))) %>%
  summarise(x = median(umap_1), y = median(umap_2), .groups = "drop")
mx <- setNames(ctr$x, ctr$node); my <- setNames(ctr$y, ctr$node)

edges_of <- function(grp) {
  conn %>%
    filter(group == grp, connectivity > 0.02) %>%
    mutate(x = mx[as.character(source)], y = my[as.character(source)],
           xend = mx[as.character(target)], yend = my[as.character(target)])
}

group_meta <- tribble(
  ~grp,    ~title,
  "joint", "Joint",
  "wt",    "Non-diabetic",
  "db",    "Diabetic"
)

paga_panel <- function(grp, title) {
  e <- edges_of(grp)
  ggplot() +
    rasterise(geom_point(aes(umap_1, umap_2), df, color = "grey90", stroke = 0, size = 0.3), dpi = 900) +
    # PAGA edges: width + opacity scale with connectivity (classic PAGA look)
    geom_segment(aes(x, y, xend = xend, yend = yend,
                     linewidth = connectivity, alpha = connectivity), e,
                 color = "grey15", lineend = "round") +
    geom_point(aes(x, y), ctr, color = "#fff", size = 1.7) +
    geom_point(aes(x, y, color = factor(node)), ctr, size = 1.2) +
    geom_shadowtext(aes(x, y, label = node), ctr, color = "#000", bg.color = "#fff",
                    size = pt2mm(fs_min - 1)) +
    scale_color_manual(values = cluster_pal, guide = "none") +
    scale_linewidth_continuous(range = c(0.1, 1.5), limits = c(0, 1), guide = "none") +
    scale_alpha_continuous(range = c(0.06, 0.95), limits = c(0, 1), guide = "none") +
    labs(title = title, x = expression("UMAP"[1]), y = expression("UMAP"[2])) +
    c(theme_embedding(arrow = arrow), list(theme(
      plot.title = element_text(size = fs_base, hjust = 0.5),
      axis.title = element_text(size = fs_min), legend.position = "none"
    )))
}

fig <- pmap(group_meta, paga_panel) %>% wrap_plots(nrow = 1)

# 124 mm wide (4.882 in), 3 panels
save_figure(output_dir, "panel_H_paga_graph", plot = fig, width = 4.882, height = 1.85)
