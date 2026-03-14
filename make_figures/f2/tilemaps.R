library(tidyverse)
library(ggbrandon)
library(patchwork)
library(ggh4x)
library(legendry)
library(colorspace)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 7))

# ----------------------------------------------------------------------
output_dir <- "figures/f2"
source("make_figures/cfg.R")
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df()

terms <- c(
  "Angiogenesis",
  "Adipogenesis",
  "Interferon Alpha Response",
  "Interferon Gamma Response",
  "Coagulation",
  "Inflammatory Response",
  "TNF-alpha Signaling via NF-kB",
  "Peroxisome",
  "Hypoxia",
  "G2-M Checkpoint",
  "Apoptosis",
  "Unfolded Protein Response"
)

res <- "results/enrichment_results/GSEA.csv" %>%
  read_csv() %>%
  mutate(
    cluster = factor(cluster, levels(df$clusterid)),
    Term = case_match(Term, "Pperoxisome" ~ "Peroxisome", .default = Term)
  ) %>%
  filter(
    `NOM p-val` < 0.05,
    # `FDR q-val` < 0.05,
    Term %in% terms
  )

# compute luminance-based text color from the fill scale
nes_lim <- max(abs(res$NES))
div_cols <- rev(diverging_hcl(256, h = c(330, 240), c = 80, l = c(35, 95)))
nes_pal <- scales::gradient_n_pal(
  div_cols,
  values = seq(0, 1, length.out = 256)
)
res <- res %>%
  mutate(
    bg_hex = nes_pal(scales::rescale(NES, from = c(-nes_lim, nes_lim))),
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

p <- res %>%
  ggplot(aes(cluster, Term)) +
  geom_tile(aes(fill = NES)) +
  geom_text(
    aes(
      label = case_when(
        `FDR q-val` < 0.001 ~ "***",
        `FDR q-val` < 0.01 ~ "**",
        `FDR q-val` < 0.05 ~ "*",
        .default = sprintf("%.2f", round(`FDR q-val`, 2))
      ),
      color = text_color
    ),
    size = 5 / .pt
  ) +
  scale_color_identity() +
  scale_fill_gradientn(
    colours = div_cols,
    values = scales::rescale(seq(-nes_lim, nes_lim, length.out = 256)),
    breaks = scales::pretty_breaks(n = 3),
    guide = guide_colorbar(
      theme = theme(
        legend.key.width  = unit(2, "lines"),
        legend.key.height = unit(0.4, "lines")
      )
    )
  ) +
  scale_y_discrete(limits = rev) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Normalized enrichment score"
  ) +
  theme(
    axis.text.y = element_text(size = 6.5),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6, hjust = 0.5),
    legend.position = "bottom",
    legend.location = "panel",
    legend.justification = "center",
    legend.title.position = "top",
    legend.text.position = "bottom",
    legend.margin = margin()
  )

save_figure(output_dir, "tilemaps", plot = p, width = 3.5, height = 1.9)
