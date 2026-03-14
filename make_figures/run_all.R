#!/usr/bin/env Rscript

# Run all figure-making scripts
# Usage: Rscript make_figures/run_all.R

scripts <- c(
  # Figure 1
  "make_figures/f1/umaps_overview.R",
  "make_figures/f1/umaps_celloracle.R",
  "make_figures/f1/heatmap_tf_activity.R",
  "make_figures/f1/stacks_timepoint.R",
  "make_figures/f1/pgd_umaps.R",

  # Figure 2
  "make_figures/f2/umaps.R",
  "make_figures/f2/barplots.R",
  "make_figures/f2/dotplots.R",
  "make_figures/f2/tilemaps.R",

  # Figure 3
  "make_figures/f3/umaps.R",
  "make_figures/f3/heatmaps_tde.R",
  "make_figures/f3/lineplots.R",

  # Figure 4
  "make_figures/f4/umaps.R",
  "make_figures/f4/heatmaps_clusters.R",
  "make_figures/f4/heatmaps_tde.R",

  # Figure 5
  "make_figures/f5/umaps_ps_markov.R",
  "make_figures/f5/heatmaps_mc_transitions.R",
  "make_figures/f5/scatterplots_ps.R",

  # Figure 6
  "make_figures/f6/umaps_ps_markov.R",
  "make_figures/f6/heatmaps_mc_transitions.R",
  "make_figures/f6/scatterplots_ps.R",

  # Supplemental figures
  "make_figures/sf1/umaps_gene_features.R",
  "make_figures/sf2/umaps_gene_features.R",
  "make_figures/sf3/umaps_ps.R",
  "make_figures/sf3/scatterplots_ps.R",
  "make_figures/sf4/barplot_go_terms.R"
)

for (s in scripts) {
  cat(sprintf(">>> %s\n", s))
  tryCatch(
    source(s, local = new.env(parent = globalenv())),
    error = function(e) cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
  )
}

cat("Done.\n")
