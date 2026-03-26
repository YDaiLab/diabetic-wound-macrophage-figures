library(tidyverse)
library(ggbrandon)
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(ggrastr)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 7))

# ----------------------------------------------------------------------
output_dir <- "figures/f4"
source("make_figures/cfg.R")
source("make_figures/load_seurat.R")

genes <- c("Brd3", "Irf9", "Hes1", "Zmynd11", "Zbtb46", "Stat5b")
class <- c("Early-stage", "", "Late-stage", "", "APC-like", "")
class <- factor(class, c("Early-stage", "Late-stage", "APC-like", ""))

class_genes <- setNames(class, genes)
class_genes_upper <- setNames(class, toupper(genes))
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
mat1 <- GetAssayData(cells, assay = "Z_avg", layer = "scale.data")

cells <- ScaleData(cells, assay = "RNA", features = genes, scale.max = 2.5)
mat2 <- GetAssayData(cells, assay = "RNA", layer = "scale.data")
# ----------------------------------------------------------------------
df <- load_df(pgd = TRUE)
dat1 <- get_data(df, mat1) %>%
  mutate(gene = factor(toupper(gene), toupper(genes)))
dat2 <- get_data(df, mat2)

p1 <- dat1 %>%
  arrange(abs(value)) %>%
  ggplot(aes(umap_1, umap_2)) +
  geom_point_rast(aes(color = value), size = 0.3, stroke = 0, raster.dpi = 600) +
  geom_shadowtext(
    aes(label = clusterid),
    df %>%
      group_by(clusterid) %>%
      summarise_if(is.numeric, median),
    color = "#000",
    bg.color = "#fff",
    size = 5 / .pt
  ) +
  facet_wrap2(
    ~ class_genes_upper[gene] + gene,
    axes = "all",
    nrow = 2,
    strip = strip_nested(
      background_x = list(
        element_rect(fill = "#EEEEEE"),
        element_rect(fill = "#EEEEEE"),
        element_rect(fill = "#EEEEEE"),
        element_blank(),
        element_rect(fill = "transparent"),
        element_rect(fill = "transparent"),
        element_rect(fill = "transparent"),
        element_rect(fill = "transparent"),
        element_rect(fill = "transparent"),
        element_rect(fill = "transparent")
      )
    ),
    dir = "h"
  ) +
  scale_color_distiller(
    palette = "RdYlBu",
    limits = c(-2.5, 2.5),
    breaks = scales::pretty_breaks(n = 3),
    guide = guide_colorbar(
      theme = theme(
        legend.key.width  = unit(0.4, "lines"),
        legend.key.height = unit(2, "lines")
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
    legend.text = element_text(size = 5),
    legend.justification = "top",
    legend.title = element_text(size = 5.5),
    legend.margin = margin(),
    strip.text = element_text(size = 6),
    ggh4x.facet.nestline = element_line(linewidth = 0.2),
    panel.spacing.y = unit(0, "lines")
  )

p2 <- dat2 %>%
  arrange(abs(value)) %>%
  ggplot(aes(umap_1, umap_2)) +
  geom_point_rast(aes(color = value), size = 0.3, stroke = 0, raster.dpi = 600) +
  geom_shadowtext(
    aes(label = clusterid),
    df %>%
      group_by(clusterid) %>%
      summarise_if(is.numeric, median),
    color = "#000",
    bg.color = "#fff",
    size = 5 / .pt
  ) +
  facet_wrap2(
    ~ class_genes[gene] + gene,
    axes = "all",
    nrow = 2,
    strip = strip_nested(
      background_x = list(
        element_rect(fill = "#EEEEEE"),
        element_rect(fill = "#EEEEEE"),
        element_rect(fill = "#EEEEEE"),
        element_blank(),
        element_rect(fill = "transparent"),
        element_rect(fill = "transparent"),
        element_rect(fill = "transparent"),
        element_rect(fill = "transparent"),
        element_rect(fill = "transparent"),
        element_rect(fill = "transparent")
      ),
      text_x = list(
        element_text(face = "plain"),
        element_text(face = "plain"),
        element_text(face = "plain"),
        element_text(face = "italic"),
        element_text(face = "italic"),
        element_text(face = "italic"),
        element_text(face = "italic"),
        element_text(face = "italic"),
        element_text(face = "italic"),
        element_text(face = "italic")
      )
    ),
    dir = "h"
  ) +
  scale_color_distiller(
    palette = "RdBu",
    limits = c(-2.5, 2.5),
    breaks = scales::pretty_breaks(n = 3),
    guide = guide_colorbar(
      theme = theme(
        legend.key.width  = unit(0.4, "lines"),
        legend.key.height = unit(2, "lines")
      )
    )
  ) +
  theme_dimred2(arrow = arrow) +
  labs(
    x = expression("UMAP"[1]),
    y = expression("UMAP"[2]),
    color = "Scaled\nexpression"
  ) +
  theme(
    axis.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5.5),
    legend.margin = margin(),
    legend.justification = "top",
    strip.text = element_text(size = 6),
    ggh4x.facet.nestline = element_line(linewidth = 0.2),
    panel.spacing.y = unit(0, "lines")
  )

# ----------------------------------------------------------------------
save_figure(output_dir, "umaps_p1", plot = p1, width = 3.25, height = 2.35)
save_figure(output_dir, "umaps_p2", plot = p2, width = 3.25, height = 2.35)
