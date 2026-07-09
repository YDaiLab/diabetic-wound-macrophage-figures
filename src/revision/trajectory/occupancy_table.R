library(tidyverse)

# Revision task T6 / R1-3 — per-cluster composition by condition x timepoint.
# Reconciles the contradictory Cluster 12 / Cluster 1 statements and quantifies
# the occupancy differences that motivate the per-condition trajectory check (T1).
#
# Lightweight post-processing of the pipeline output, so it lives here in moma
# (not in the macrophages workflow). Reads the cells metadata produced by
# src/4_1-extract_cells_meta.R.
# ----------------------------------------------------------------------
input_file <- "results/cells_meta.csv"
out_long <- "results/revision/occupancy_long.csv"
out_wide <- "results/revision/occupancy_wide.csv"
out_origin <- "results/revision/occupancy_origin.csv"
# ----------------------------------------------------------------------
dir.create("results/revision", recursive = TRUE, showWarnings = FALSE)

meta <- read_csv(input_file) %>%
  transmute(
    cell = Row.names,
    clusterid = fct_inseq(as.character(clusterid)),
    condition = factor(condition, c("Non-diabetic", "Diabetic")),
    timepoint = factor(timepoint, c("D3", "D6", "D10"))
  )

long <- meta %>%
  count(clusterid, condition, timepoint, name = "n") %>%
  complete(clusterid, condition, timepoint, fill = list(n = 0)) %>%
  group_by(clusterid) %>%
  mutate(frac_of_cluster = n / sum(n)) %>%
  ungroup() %>%
  group_by(condition, timepoint) %>%
  mutate(frac_of_condition_timepoint = n / sum(n)) %>%
  ungroup() %>%
  arrange(clusterid, condition, timepoint)

write_csv(long, out_long)

# What each cluster is made of (answers R1-3): rows = clusters, columns =
# condition_timepoint, values = fraction of that cluster.
wide <- long %>%
  mutate(group = paste(condition, timepoint)) %>%
  select(clusterid, group, frac_of_cluster) %>%
  pivot_wider(names_from = group, values_from = frac_of_cluster) %>%
  left_join(count(meta, clusterid, name = "n_total"), by = "clusterid") %>%
  arrange(clusterid)

write_csv(wide, out_wide)

# Origin (D3) cells per condition — the T1 stability check (infer_tree uses
# origin.celltype = "D3"; a condition thin at D3 gives an unstable tree origin).
origin <- meta %>%
  count(condition, timepoint, name = "n_cells") %>%
  complete(condition, timepoint, fill = list(n_cells = 0)) %>%
  pivot_wider(names_from = timepoint, values_from = n_cells)

write_csv(origin, out_origin)

message("Cells per condition x timepoint:")
print(origin)
