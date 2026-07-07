library(tidyverse)
source("assets/theme_nature.R")
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(colorspace)
library(ggnewscale)
library(Seurat)

theme_set(theme_nature())

# ----------------------------------------------------------------------
output_dir <- "figures/f2"
source("make_figures/cfg.R")
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df()

cells <- readRDS("box/results/seurat/cells.rds")
de_results <- "results/de_tests.xlsx" %>%
  readxl::read_excel(sheet = "one_vs_rest") %>%
  mutate(
    cluster = factor(cluster, levels(df$clusterid))
  )

genes <- de_results %>%
  filter(p_val_adj < 0.05, pct.1 > 0.3) %>%
  slice_max(order_by = avg_log2FC, by = cluster, n = 5) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  pull(gene) %>%
  unique()

dotplot <- DotPlot(cells, features = genes)
dat <- as_tibble(dotplot@data) %>%
  mutate(id = factor(id, levels(df$clusterid))) %>%
  left_join(
    de_results %>%
      filter(avg_log2FC > 0, pct.1 > 0.3) %>%
      select(gene, cluster, p_val_adj) %>%
      mutate(sig = case_when(
        p_val_adj < 0.001 ~ "***",
        p_val_adj < 0.01 ~ "**",
        p_val_adj < 0.05 ~ "*",
        .default = ""
      )),
    by = c("features.plot" = "gene", "id" = "cluster")
  ) %>%
  mutate(features.plot = factor(features.plot, levels = genes))

# map avg.exp.scaled to the RdBu palette and compute luminance
lim <- max(abs(dat$avg.exp.scaled))
rdbu_pal <- scales::gradient_n_pal(
  rev(RColorBrewer::brewer.pal(11, "RdBu")),
  values = scales::rescale(seq(-lim, lim, length.out = 11))
)
dat <- dat %>%
  mutate(
    bg_hex = rdbu_pal(scales::rescale(avg.exp.scaled, from = c(-lim, lim))),
    text_color = ifelse(
      colorspace::hex2RGB(bg_hex) %>%
        as("polarLUV") %>%
        slot("coords") %>%
        {
          .[, "L"]
        } < 50,
      "white", "black"
    )
  )

p <- dat %>%
  ggplot(aes(id, features.plot)) +
  geom_point(aes(color = avg.exp.scaled, size = pct.exp), stroke = 0) +
  scale_color_distiller(
    palette = "RdBu",
    limits = c(-1, 1) * lim,
    breaks = scales::pretty_breaks(n = 3),
    guide = guide_colorbar(
      title = "Average expression scaled",
      order = 2,
      theme = theme(
        legend.key.width  = unit(2, "lines"),
        legend.key.height = unit(0.4, "lines")
      )
    )
  ) +
  # new_scale_color() +
  # geom_text(
  #   aes(
  #     label = sig,
  #     color = text_color
  #   ),
  #   vjust = 0.75, size = pt2mm(fs_min)
  # ) +
  # scale_color_identity() +
  scale_size_area(
    labels = scales::percent_format(scale = 1),
    max_size = 4,
    breaks = seq(25, 75, by = 25),
    guide = guide_legend(order = 1)
  ) +
  scale_y_discrete(limits = rev) +
  theme(
    axis.text.y = element_text(face = "italic", size = fs_small),
    legend.text = element_text(size = fs_small, hjust = 0.5),
    legend.title = element_text(size = fs_small, hjust = 0.5),
    legend.position = "bottom",
    legend.location = "panel",
    legend.justification = "left",
    legend.title.position = "top",
    legend.text.position = "bottom",
    legend.margin = margin(),
    legend.key.height = unit(0, "lines")
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Average expression scaled",
    size = "Percent expressing"
  )

save_figure(output_dir, "dotplots", plot = p, width = 3.3, height = 6)
