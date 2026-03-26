source("make_tables/cfg.R")

df_500 <-
  "box/results/pgd_celloracle_overexpress/mc_transitions/500.parquet" %>%
  arrow::read_parquet()

df_0 <-
  "box/results/pgd_celloracle_overexpress/mc_transitions/0.parquet" %>%
  arrow::read_parquet()

mat <- bind_rows(
  df_0 %>%
    mutate(step = 0),
  df_500 %>%
    mutate(step = 500)
) %>%
  count(step, factor, clusterid) %>%
  pivot_wider(names_from = step, values_from = n) %>%
  mutate(
    delta = `500` - `0`,
    pct_change = delta / `0` * 100
  ) %>%
  select(factor, clusterid, pct_change)

output_file <- get_table_path("Table S8", "markov_overexpression.csv")
write_csv(mat, output_file)
