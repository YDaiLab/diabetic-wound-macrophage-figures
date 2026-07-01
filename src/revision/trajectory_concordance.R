library(tidyverse)

# Revision task T1 / R2-2 (editor-flagged) — quantify agreement between the joint
# trajectory and the independently-inferred per-condition (wt, db) trajectories.
# This is the number backing the per-condition figure: a visual "the trajectory
# is still there" is not enough for the editor-flagged point.
#
#   1. Pseudotime concordance — Spearman(per-condition pt, joint pt) on the
#      condition's cells, overall and within each joint branch.
#   2. Branch crosstab       — do cells forming a joint branch also form a
#      coherent branch in the per-condition tree? (occupancy differs, topology
#      shouldn't).
#   3. Cluster ordering      — sequence of joint cluster labels along each
#      per-condition branch (topology check robust to per-condition k-means
#      using different node numbers).
#
# Per-condition k-means cluster numbers are NOT assumed to match the joint
# numbering; everything is anchored on cell ids and the joint clusterid.
# Lightweight, so it lives in moma. Reads the .rds trees synced from macrophages.
# ----------------------------------------------------------------------
joint_file <- "box/results/lamian/infer_tree.rds"
cond_files <- c(
  wt = "box/results/revision/lamian/wt/infer_tree.rds",
  db = "box/results/revision/lamian/db/infer_tree.rds"
)
out_pt <- "results/revision/concordance_pseudotime.csv"
out_branch <- "results/revision/concordance_branch_crosstab.csv"
out_order <- "results/revision/concordance_cluster_order.csv"
# ----------------------------------------------------------------------
dir.create("results/revision", recursive = TRUE, showWarnings = FALSE)

# Per-cell branch membership from a Lamian `order` list (branch -> ordered cells).
cell_branch <- function(res) {
  imap(res[["order"]], ~ tibble(cell = .x, branch = .y)) %>%
    bind_rows()
}

# Lamian's `pseudotime`/`clusterid` list shared backbone cells once per branch,
# so names repeat. Dedupe to one row per cell for the pseudotime/cluster tables
# (branch membership is intentionally kept multi-valued via cell_branch()).
uniq_cell <- function(x) distinct(x, cell, .keep_all = TRUE)

joint <- readRDS(joint_file)
joint_pt <- enframe(joint[["pseudotime"]], name = "cell", value = "pt_joint") %>%
  uniq_cell()
joint_cl <- enframe(joint[["clusterid"]], name = "cell", value = "cluster_joint") %>%
  uniq_cell() %>%
  mutate(cluster_joint = fct_inseq(as.character(cluster_joint)))
joint_branch <- cell_branch(joint) %>% rename(branch_joint = branch)

pt_rows <- list()
branch_rows <- list()
order_rows <- list()

for (cond in names(cond_files)) {
  res <- readRDS(cond_files[[cond]])
  cond_pt <- enframe(res[["pseudotime"]], name = "cell", value = "pt_cond") %>%
    uniq_cell()
  cond_branch <- cell_branch(res) %>% rename(branch_cond = branch)

  # One row per shared cell: pseudotime concordance is a per-cell comparison.
  pt_tab <- cond_pt %>%
    inner_join(joint_pt, by = "cell") %>%
    inner_join(joint_cl, by = "cell")

  overall <- tibble(
    condition = cond, branch_joint = "ALL", n = nrow(pt_tab),
    spearman = cor(pt_tab$pt_cond, pt_tab$pt_joint, method = "spearman")
  )
  per_branch <- pt_tab %>%
    inner_join(joint_branch, by = "cell") %>%
    group_by(branch_joint) %>%
    filter(n() >= 10) %>%
    summarise(
      condition = cond, n = n(),
      spearman = cor(pt_cond, pt_joint, method = "spearman"),
      .groups = "drop"
    )
  pt_rows[[cond]] <- bind_rows(overall, per_branch) %>% relocate(condition)

  # Branch crosstab: joint-branch membership x per-condition-branch membership.
  branch_rows[[cond]] <- joint_branch %>%
    inner_join(cond_branch, by = "cell") %>%
    count(condition = cond, branch_joint, branch_cond, name = "n")

  order_rows[[cond]] <- pt_tab %>%
    inner_join(cond_branch, by = "cell") %>%
    group_by(branch_cond) %>%
    arrange(pt_cond, .by_group = TRUE) %>%
    mutate(condition = cond, pt_decile = ntile(pt_cond, 10)) %>%
    group_by(condition, branch_cond, pt_decile) %>%
    summarise(
      n = n(),
      modal_joint_cluster = names(which.max(table(cluster_joint))),
      .groups = "drop"
    )
}

write_csv(bind_rows(pt_rows), out_pt)
write_csv(bind_rows(branch_rows), out_branch)
write_csv(bind_rows(order_rows), out_order)

message("Pseudotime concordance (Spearman vs joint):")
print(bind_rows(pt_rows) %>% filter(branch_joint == "ALL"))
