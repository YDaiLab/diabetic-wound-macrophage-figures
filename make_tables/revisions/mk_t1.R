source("make_tables/cfg.R")

# Revision R1-3 / R1-5 — cluster composition by genotype and time point.
#
# A page-formatted display item (12 rows x 8 columns), not a spreadsheet: the
# point is that a reader can check the cluster descriptions in the Results while
# reading them, which is what the Reviewer found unclear.
#
# Percentages are of the cluster, the denominator those statements use, and the
# six columns sum to 100 per row, so the marginals the Results quote are read by
# adding the relevant cells: cluster 12 is 68.8 + 20.3 = 89% day 3 and
# 68.8 + 7.3 + 0.3 = 76% non-diabetic; cluster 1 is 34.8 + 25.8 = 61% diabetic
# day 6 and day 10.
#
# Reads the lightweight summary written by src/revision/trajectory/occupancy_table.R.
# ----------------------------------------------------------------------
input_file <- "results/revision/occupancy_long.csv"
group_levels <- c(
  "Non-diabetic D3", "Non-diabetic D6", "Non-diabetic D10",
  "Diabetic D3", "Diabetic D6", "Diabetic D10"
)
# ----------------------------------------------------------------------

long <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(group = factor(paste(condition, timepoint), group_levels))

tbl <- long %>%
  transmute(clusterid, group, pct = round(100 * frac_of_cluster, 1)) %>%
  pivot_wider(names_from = group, values_from = pct) %>%
  left_join(
    long %>% group_by(clusterid) %>% summarise(n = sum(n), .groups = "drop"),
    by = "clusterid"
  ) %>%
  rename(Cluster = clusterid, `Cells in cluster` = n) %>%
  arrange(Cluster)

# the six percentages must account for the whole cluster
stopifnot(
  nrow(tbl) == 12,
  all(abs(rowSums(tbl[, group_levels]) - 100) < 0.2),
  sum(tbl$`Cells in cluster`) == 6109
)

output_file <- get_table_path(1, "cluster_composition.csv", kind = "Table")
write_csv(tbl, output_file)
message("wrote ", output_file)
print(tbl, width = Inf)
