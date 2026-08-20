library(tidyverse)
library(arrow)

# Revision R2-5, panel d — do the biological conclusions depend on PGD?
#
# CellOracle was run twice on the same cells, the same GRN and the same
# perturbations, differing only in the embedding the vector field is projected
# onto: the conventional UMAP (results/celloracle) and the PGD layout
# (results/pgd_celloracle). This aggregates both to the quantity the manuscript
# actually reports in Fig. 5a,b — summed negative (inhibitory) and positive
# (promoting) perturbation scores per TR per trajectory branch.
#
# NB compare RANKS, not magnitudes. CellOracle fits the vector-field grid to
# each embedding separately (8,598 vs 5,013 grid points per TR), so the sums are
# on different scales by construction.
# ----------------------------------------------------------------------
# Both perturbation directions the manuscript reports: knockout (Fig. 5) and
# overexpression (Fig. 6). Panel c plots the knockout half; the overexpression
# concordance is quoted in prose, so it lives in the same table rather than a
# second figure.
runs <- tribble(
  ~perturbation, ~umap_file, ~pgd_file,
  "Knockout",
  "box/results/celloracle/perturb_scores.parquet",
  "box/results/pgd_celloracle/perturb_scores.parquet",
  "Overexpression",
  "box/results/celloracle_overexpress/perturb_scores.parquet",
  "box/results/pgd_celloracle_overexpress/perturb_scores.parquet"
)
output_file <- "results/revision/pgd_perturb_comparison.csv"
# ----------------------------------------------------------------------
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Grid cells outside any branch carry no branch-level score. is.na(subbranch)
# keeps the branch-level rows only, as make_figures/f5/scatterplots_ps.R does:
# the subbranch rows re-partition the SAME grid points into segments (12_6, 6_7,
# ...), and junction points fall in more than one segment, so summing both levels
# double-counts unevenly.
branch_sums <- function(fname) {
  read_parquet(fname) %>%
    filter(!is.na(branch), is.na(subbranch)) %>%
    group_by(factor, branch) %>%
    summarise(
      neg = sum(pmin(score, 0)),
      pos = sum(pmax(score, 0)),
      .groups = "drop"
    )
}

cmp <- runs %>%
  pmap_dfr(function(perturbation, umap_file, pgd_file) {
    inner_join(
      branch_sums(umap_file), branch_sums(pgd_file),
      by = c("factor", "branch"), suffix = c("_umap", "_pgd")
    ) %>%
      mutate(perturbation = perturbation)
  }) %>%
  group_by(perturbation, branch) %>%
  mutate(
    # rank 1 = strongest inhibition / strongest promotion within a branch
    rank_neg_umap = rank(neg_umap), rank_neg_pgd = rank(neg_pgd),
    rank_pos_umap = rank(-pos_umap), rank_pos_pgd = rank(-pos_pgd)
  ) %>%
  ungroup() %>%
  select(perturbation, everything()) %>%
  arrange(perturbation, branch, rank_neg_umap)

stopifnot(
  nrow(cmp) ==
    n_distinct(cmp$factor) * n_distinct(cmp$branch) * n_distinct(cmp$perturbation)
)
write_csv(cmp, output_file)

# Per-branch concordance — the numbers quoted in the response letter.
concordance <- cmp %>%
  group_by(perturbation, branch) %>%
  summarise(
    rho_neg = cor(neg_umap, neg_pgd, method = "spearman"),
    rho_pos = cor(pos_umap, pos_pgd, method = "spearman"),
    top10_shared = length(intersect(
      factor[rank_neg_umap <= 10], factor[rank_neg_pgd <= 10]
    )),
    .groups = "drop"
  )

interpreted <- c("branch: 12,6,7,5", "branch: 12,6,9,3,10", "branch: 12,6,9,3,4,8")
concordance %>%
  filter(branch %in% interpreted) %>%
  group_by(perturbation) %>%
  summarise(
    inhibitory = sprintf("%.2f-%.2f", min(rho_neg), max(rho_neg)),
    promoting = sprintf("%.2f-%.2f", min(rho_pos), max(rho_pos)),
    .groups = "drop"
  ) %>%
  print()
print(concordance, n = Inf)
