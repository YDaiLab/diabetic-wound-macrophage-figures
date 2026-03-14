library(tidyverse)

# ----------------------------------------------------------------------
assets_dir <- "assets"
# ----------------------------------------------------------------------
load_pal <- function(fname) {
  require(dplyr)
  require(readr)
  require(tibble)

  readr::read_csv(fname) %>%
    dplyr::mutate(`...1` = as.character(`...1`)) %>% # nolint
    tibble::deframe()
}

gen_pal <- function(df, group_col, seed = 42) {
  require(dplyr)
  require(Palo)
  require(scales)

  set.seed(seed)
  Palo::Palo(
    as.matrix(select(df, starts_with("umap"))),
    df[[group_col]],
    scales::hue_pal()(n_distinct(df[[group_col]]))
  )
}

load_df <- function(pgd = FALSE) {
  require(forcats)

  fname <- if (pgd) {
    "results/pgd_cells_meta.csv"
  } else {
    "results/cells_meta.csv"
  }

  df <- fname %>%
    read_csv() %>%
    mutate(
      sample = fct(sample, names(sample_pal)),
      clusterid = fct_inseq(as.character(clusterid)) # nolint
    ) %>%
    sample_frac(1L)

  cluster_order <- df %>%
    count(clusterid, condition) %>% # nolint
    group_by(clusterid) %>%
    mutate(p = n / sum(n)) %>%
    filter(condition == "Diabetic") %>%
    arrange(p) %>% # nolint
    pull(clusterid)

  df %>%
    mutate(
      clusterid = factor(clusterid, cluster_order), # nolint
      condition = factor(condition, c("Non-diabetic", "Diabetic")) # nolint
    )
}
# ----------------------------------------------------------------------
dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)

sample_pal <- c(
  "Non-diabetic D3" = "#A7DBD8",
  "Non-diabetic D6" = "#1E988A",
  "Non-diabetic D10" = "#004D44",
  "Diabetic D3" = "#FDDCAA",
  "Diabetic D6" = "#F5A623",
  "Diabetic D10" = "#D4700A"
)

condition_pal <- c(
  "Non-diabetic" = "#1E988A",
  "Diabetic" = "#F5A623"
)

# cluster palette
fname <- file.path(assets_dir, "cluster_pal.csv")
if (file.exists(fname)) {
  cluster_pal <- load_pal(fname)
} else {
  cluster_pal <- gen_pal(load_df(), "clusterid", seed = 42)
  write.csv(cluster_pal, fname)
}

arrow <-
  grid::arrow(type = "closed", angle = 25, length = grid::unit(0.1, "lines"))

save_figure <- function(output_dir, filename, plot = last_plot(),
                        width, height, dpi = 900) {
  for (fmt in c("png", "pdf", "svg")) {
    fmt_dir <- file.path(output_dir, fmt)
    dir.create(fmt_dir, recursive = TRUE, showWarnings = FALSE)
    out <- file.path(fmt_dir, paste0(filename, ".", fmt))
    if (fmt == "png") {
      ggsave(out, plot, width = width, height = height, dpi = dpi)
    } else {
      ggsave(out, plot, width = width, height = height)
    }
  }
}
