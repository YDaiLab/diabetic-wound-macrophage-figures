# Revision R2-2 — export branch orders + full PCA for every raw draw, so PGD can
# be run per draw (src/revision/run_pgd_draws.py).
#   box/results/revision/lamian_condpca_seeds_raw/{wt,db}/draw{1..10}/infer_tree.rds
# ->  results/revision/pgd_inputs_draws/{cond}_draw{NN}_{orders,pca}.csv
library(tidyverse)

draw_dir <- "box/results/revision/lamian_condpca_draws"
outdir <- "results/revision/pgd_inputs_draws"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (cond in c("wt", "db")) {
  for (d in 1:10) {
    res <- readRDS(file.path(draw_dir, cond, paste0("seed", d), "infer_tree.rds"))
    orders <- imap_dfr(res$order, ~ tibble(branch = .y, pos = seq_along(.x), cell = .x))
    write_csv(orders, file.path(outdir, sprintf("%s_draw%02d_orders.csv", cond, d)))
    pca <- as.data.frame(res$pca_full) %>% rownames_to_column("cell")
    write_csv(pca, file.path(outdir, sprintf("%s_draw%02d_pca.csv", cond, d)))
  }
  message(cond, ": exported 10 draws")
}
