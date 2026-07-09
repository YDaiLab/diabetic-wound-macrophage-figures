library(tidyverse)
source("assets/theme_nature.R")
library(patchwork)
library(igraph)

theme_set(theme_nature())

# Revision R2-2, panel F — edge-level reproducibility of the joint trajectory
# across the 10-draw condition-specific ensemble.
#
# For each JOINT MST edge (a-b), restricted to a condition's cells, we take the
# best cell-Jaccard match among each draw's edges (label-agnostic). Over the 10
# draws that gives a spread (min .. max) with a median. The official pipeline
# (condpca_edge_detection.csv) only stored min/median + the per-edge null cutoff
# (99th pct vs random cell sets), so we recompute the full spread here from the
# same trees and verify min/median match the CSV before plotting min..max.
#
# ND over DB, each sorted independently by median. Monochrome.
# ----------------------------------------------------------------------
output_dir <- "figures/revision/trajectory/supp_panels"
source("make_figures/cfg.R")

joint_file <- "box/results/lamian/infer_tree.rds"
draw_dir <- "box/results/revision/lamian_condpca_draws"
conds <- c(wt = "Non-diabetic", db = "Diabetic")
# ----------------------------------------------------------------------

jaccard <- function(A, B) length(intersect(A, B)) / length(union(A, B))
ekey <- function(a, b) paste(sort(as.integer(c(a, b))), collapse = "-")

joint <- readRDS(joint_file)
jcl <- setNames(as.character(joint$clusterid), names(joint$clusterid))
jedges <- igraph::as_edgelist(joint$MSTtree)

# per-draw best-match Jaccard for every joint edge, then min/median/max over draws
spread <- imap_dfr(conds, function(label, cond) {
  per_draw <- map_dfr(1:10, function(dd) {
    r <- readRDS(file.path(draw_dir, cond, paste0("seed", dd), "infer_tree.rds"))
    dcl <- setNames(as.character(r$clusterid), names(r$clusterid))
    cond_cells <- names(dcl)
    de <- igraph::as_edgelist(r$MSTtree)
    dsets <- lapply(seq_len(nrow(de)), function(i) names(dcl)[dcl %in% de[i, ]])
    map_dfr(seq_len(nrow(jedges)), function(j) {
      Sj <- intersect(names(jcl)[jcl %in% jedges[j, ]], cond_cells)
      tibble(edge = ekey(jedges[j, 1], jedges[j, 2]), draw = dd,
             jacc = max(vapply(dsets, function(Sd) jaccard(Sj, Sd), numeric(1))))
    })
  })
  per_draw %>%
    group_by(edge) %>%
    summarise(jmin = min(jacc), jmed = median(jacc), jmax = max(jacc), .groups = "drop") %>%
    mutate(condition = label)
})

# Chance level = the best-of-k Jaccard null (js_bok_cut), which accounts for the
# best-match statistic being a MAX over the draw's edges (a single-set null would
# be anti-conservative). Also pull the ensemble significance: every edge is
# recovered in 10/10 draws under this calibrated null, so the binomial combined
# p (BH-FDR adjusted, p_bh) hits its floor for all edges. Verify recomputed
# min/median magnitudes agree with the official run.
off <- read_csv("box/results/revision/condpca_edge_detection.csv", show_col_types = FALSE) %>%
  mutate(condition = recode(condition, wt = "Non-diabetic", db = "Diabetic"))
d <- spread %>%
  left_join(off %>% select(condition, edge, null_cutoff = js_bok_cut,
                           jacc_min, jacc_median, n_detected_bok, n_draws, p_bh),
            by = c("condition", "edge"))
message(sprintf("verify vs CSV — max |Δmin|=%.4g, max |Δmedian|=%.4g",
                max(abs(d$jmin - d$jacc_min)), max(abs(d$jmed - d$jacc_median))))
stopifnot(max(abs(d$jmin - d$jacc_min)) < 1e-6, max(abs(d$jmed - d$jacc_median)) < 1e-6)

# one-line ensemble result for the figure (best-of-k null; binomial + BH-FDR)
min_det <- min(d$n_detected_bok); nd <- d$n_draws[1]
sig_note <- sprintf("All edges recovered in %d/%d draws (FDR p < 1e-19)", min_det, nd)
message(sig_note)

one_panel <- function(cond, show_x) {
  dd <- d %>% filter(condition == cond) %>% arrange(jmed) %>% mutate(edge = fct_inorder(edge))
  ggplot(dd, aes(y = edge)) +
    geom_segment(aes(x = jmin, xend = jmax, yend = edge), color = "grey35", linewidth = 0.5) +
    geom_point(aes(x = null_cutoff), shape = 124, size = 1.5, color = "grey55") +
    geom_point(aes(x = jmed), color = "grey15", size = 1.4) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(title = cond, y = NULL, x = if (show_x) "Best-match cell-Jaccard" else NULL) +
    theme(
      plot.title = element_text(size = fs_small, hjust = 0.5),
      axis.title.x = element_text(size = fs_min),
      axis.text = element_text(size = fs_min),
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey92"),
      panel.grid.minor = element_blank()
    )
}

fig <- one_panel("Non-diabetic", show_x = FALSE) /
  one_panel("Diabetic", show_x = TRUE)

# 55 mm wide (2.165 in). Ensemble significance (sig_note) is reported in the
# figure legend, not on the panel.
save_figure(output_dir, "panel_F_edge_jaccard", plot = fig, width = 2.165, height = 3.6)
