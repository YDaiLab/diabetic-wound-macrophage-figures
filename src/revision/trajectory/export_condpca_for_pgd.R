# Revision R2-2 — export what the per-condition PGD run (Python) needs from the
# variant-C (condition-specific PCA) trajectories: each condition's branch orders
# and its OWN condition-specific PC embedding. Written as plain CSVs so the Python
# side needs no rpy2.
#
#   results/revision/pgd_inputs/{wt,db}_orders.csv  (branch, pos, cell)
#   results/revision/pgd_inputs/{wt,db}_pca.csv     (cell, PC1..PCk)
library(tidyverse)

cond_files <- c(
  wt = "box/results/revision/lamian_condpca/wt/infer_tree.rds",
  db = "box/results/revision/lamian_condpca/db/infer_tree.rds"
)
outdir <- "results/revision/pgd_inputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (cond in names(cond_files)) {
  res <- readRDS(cond_files[[cond]])

  # Branch orders: values of res$order are the ordered cell ids per branch.
  # (Lamian repeats shared backbone cells across branches; build_graph handles it.)
  orders <- imap_dfr(res$order, ~ tibble(branch = .y, pos = seq_along(.x), cell = .x))
  write_csv(orders, file.path(outdir, paste0(cond, "_orders.csv")))

  # Full (npcs-dim) condition-specific PC embedding, NOT the elbow-truncated
  # res$pca that infer_tree_structure keeps for the tree. PGD diffuses on the
  # full embedding so the layout is not constrained to the trajectory's dims.
  pca <- as.data.frame(res$pca_full) %>% rownames_to_column("cell")
  write_csv(pca, file.path(outdir, paste0(cond, "_pca.csv")))

  message(sprintf("[%s] %d cells, %d PCs (pca_full), %d branches",
                  cond, nrow(pca), ncol(pca) - 1, length(res$order)))
}
