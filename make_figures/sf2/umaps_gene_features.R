library(tidyverse)
library(ggbrandon)
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(ggrastr)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 7))

# ----------------------------------------------------------------------
output_dir <- "figures/sf2"
source("make_figures/cfg.R")
source("make_figures/load_seurat.R")

genes <- c(
  "Nfe2l2",
  "Phf10",
  "Nr1h3",
  "Ezh2",
  "Brd2",
  "Nr3c1"
)
class <- c(
  "Cluster 10",
  "Cluster 10",
  "Cluster 1",
  "Cluster 1",
  "Cluster 4/7/8",
  "Cluster 4/7/8"
)
class <- fct_inorder(class)

class_genes <- setNames(class, genes)
# ----------------------------------------------------------------------
get_data <- function(df, mat) {
  require(dplyr)
  require(tibble)
  require(tidyr)

  mat %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(where(is.numeric)) %>%
    left_join(
      df %>%
        select(Row.names, umap_1, umap_2, clusterid), # nolint
      by = join_by(name == Row.names) # nolint
    ) %>%
    mutate(gene = factor(gene, genes)) # nolint
}
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
cells <- load_cells()
cells <- ScaleData(cells, assay = "Z_avg", features = genes, scale.max = 2.5)
mat <- GetAssayData(cells, assay = "Z_avg", layer = "scale.data")
# ----------------------------------------------------------------------
df <- load_df()
df_pgd <- load_df(pgd = T)
dat <- get_data(df, mat)
dat_pgd <- get_data(df_pgd, mat)

p1 <- dat %>%
  arrange(abs(value)) %>%
  ggplot(aes(umap_1, umap_2)) +
  geom_point_rast(aes(color = value), size = 0.5, stroke = 0, raster.dpi = 600) +
  geom_shadowtext(
    aes(label = clusterid),
    df %>%
      group_by(clusterid) %>%
      summarise_if(is.numeric, median),
    color = "#000",
    bg.color = "#fff",
    size = 5 / .pt
  ) +
  facet_nested_wrap(~ class_genes[gene] + gene,
    axes = "all",
    nrow = 1,
    strip = strip_nested(
      text_x = list(
        element_text(face = "plain"),
        element_text(face = "plain")
      ),
      by_layer_x = TRUE
    )
  ) +
  scale_color_distiller(
    palette = "RdYlBu",
    limits = c(-2.5, 2.5),
    breaks = scales::pretty_breaks(n = 3),
    guide = guide_colorbar(
      theme = theme(
        legend.key.height = unit(1.5, "lines"),
        legend.key.width = unit(0.4, "lines")
      )
    )
  ) +
  theme_dimred2(arrow = arrow) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    color = "Scaled\nactivity"
  ) +
  theme(
    axis.title = element_text(size = 5),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.margin = margin(),
    strip.text = element_text(size = 7),
    ggh4x.facet.nestline = element_line(linewidth = 0.2)
  )

p2 <- dat_pgd %>%
  arrange(abs(value)) %>%
  ggplot(aes(umap_1, umap_2)) +
  geom_point_rast(aes(color = value), size = 0.5, stroke = 0, raster.dpi = 600) +
  geom_shadowtext(
    aes(label = clusterid),
    df_pgd %>%
      group_by(clusterid) %>%
      summarise_if(is.numeric, median),
    color = "#000",
    bg.color = "#fff",
    size = 5 / .pt
  ) +
  facet_nested_wrap(~ class_genes[gene] + gene,
    axes = "all",
    nrow = 1,
    strip = strip_nested(
      text_x = list(
        element_text(face = "plain"),
        element_text(face = "plain")
      ),
      by_layer_x = TRUE
    )
  ) +
  scale_color_distiller(
    palette = "RdYlBu",
    limits = c(-2.5, 2.5),
    breaks = scales::pretty_breaks(n = 3),
    guide = guide_colorbar(
      theme = theme(
        legend.key.height = unit(1.5, "lines"),
        legend.key.width = unit(0.4, "lines")
      )
    )
  ) +
  theme_dimred2(arrow = arrow) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    color = "Scaled\nactivity"
  ) +
  theme(
    axis.title = element_text(size = 5),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.margin = margin(),
    strip.text = element_text(size = 7),
    ggh4x.facet.nestline = element_line(linewidth = 0.2)
  )

combined <- p1 + plot_spacer() + p2 +
  plot_layout(guides = "collect", ncol = 1, heights = c(1, 0.1, 1)) &
  theme(
    plot.margin = margin(),
    legend.justification = "top",
    panel.spacing.y = unit(1, "lines")
  )

save_figure(output_dir, "umaps_gene_features", plot = combined, width = 7.1, height = 3.5)
