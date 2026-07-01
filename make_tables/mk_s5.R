source("make_tables/cfg.R")
library(Lamian)

res <- readRDS("box/results/lamian/z/tde.rds")

df <- lapply(names(res), function(x) {
  tmp <-
    getPopulationFit(res[[x]], type = "time")
  mat <- t(scale(t(tmp)))
  peak_time <- apply(mat, 1, which.max)

  res[[x]][["statistics"]] |>
    as.data.frame() |>
    rownames_to_column("gene") |>
    mutate(branch_name = x) |>
    mutate(
      peak_time = peak_time[gene],
      peak_time_group = case_when(
        peak_time < 10 ~ "Early",
        peak_time > 990 ~ "Late",
        .default = "Middle"
      )
    )
}) |>
  bind_rows() |>
  as_tibble()

output_file <- get_table_path(5, "trajectory_da_trs.csv")
write_csv(df, output_file)
