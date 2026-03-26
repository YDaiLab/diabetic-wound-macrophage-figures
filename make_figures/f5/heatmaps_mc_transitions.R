library(tidyverse)
library(ggbrandon)
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(ggrepel)
library(colorspace)
library(ComplexHeatmap)
library(circlize)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 7))
ht_opt$annotation_use_raster <- TRUE
options(ggrepel.max.overlaps = 6)
# ----------------------------------------------------------------------
output_dir <- "figures/f5"
source("make_figures/cfg.R")
rep_factors <- c("Cebpa", "Nr1h3", "Irf4", "Kdm1a")
rep_factors_pal <- setNames(
  RColorBrewer::brewer.pal(length(rep_factors), "Dark2"),
  rep_factors
)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df()
df_500 <-
  "box/results/pgd_celloracle/mc_transitions/500.parquet" %>%
  arrow::read_parquet()

df_0 <-
  "box/results/pgd_celloracle/mc_transitions/0.parquet" %>%
  arrow::read_parquet()

mat <- bind_rows(
  df_0 %>%
    mutate(step = 0),
  df_500 %>%
    mutate(step = 500)
) %>%
  count(step, factor, clusterid) %>%
  pivot_wider(names_from = step, values_from = n, values_fill = 0) %>%
  mutate(
    delta = `500` - `0`,
    pct_change = delta / `0` * 100
  ) %>%
  select(factor, clusterid, pct_change) %>%
  pivot_wider(names_from = clusterid, values_from = pct_change) %>%
  column_to_rownames("factor") %>%
  as.matrix()

colors <- case_match(
  rownames(mat),
  rep_factors ~ rep_factors_pal[rownames(mat)],
  .default = "black"
)
faces <- case_match(
  rownames(mat),
  rep_factors ~ "bold.italic",
  .default = "italic"
)

ht_opt$COLOR <- pals::coolwarm(100)

ht <- Heatmap(
  t(mat),
  name = "ht",
  row_order = levels(df$clusterid),
  row_names_centered = TRUE,
  row_names_rot = 0,
  row_names_gp = gpar(fontsize = 6.5),
  column_names_gp = gpar(col = colors, fontsize = 5.75, fontface = faces),
  column_dend_gp = gpar(lwd = 0.6),
  column_dend_height = unit(0.5, "lines"),
  heatmap_legend_param = list(
    title = "Percent change cluster counts\n",
    grid_height = unit(0.6, "lines"),
    legend_width = unit(4, "lines"),
    direction = "horizontal",
    at = c(-400, 0, 400),
    labels = paste0("\n", c(-400, 0, 400), "%"),
    title_gp = gpar(fontsize = 6, lineheight = 0.45),
    labels_gp = gpar(fontsize = 6),
    title_position = "lefttop"
  ),
  column_names_side = "bottom",
  column_dend_side = "top",
  row_names_side = "left",
)

hmap <-
  grid.grabExpr({
    ht_drawn <- draw(ht,
      background = "transparent",
      raster_quality = 5,
      merge_legends = TRUE,
      heatmap_legend_side = "bottom"
    )
    col_ord <- column_order(ht_drawn)
    col_names_ordered <- colnames(t(mat))[col_ord]
    rep_idx <- which(col_names_ordered %in% rep_factors)
    decorate_heatmap_body("ht", {
      nc <- ncol(t(mat))
      groups <- split(rep_idx, cumsum(c(1, diff(rep_idx) != 1)))
      for (grp in groups) {
        grid.rect(
          x = unit((min(grp) - 1) / nc, "npc"),
          y = unit(0, "npc"),
          width = unit(length(grp) / nc, "npc"),
          height = unit(1, "npc"),
          just = c("left", "bottom"),
          gp = gpar(col = "black", lwd = 1, fill = NA)
        )
      }
    })
  })

save_figure(output_dir, "heatmaps_mc_transitions", plot = hmap, width = 7.2, height = 2.4)
