# Embedding / dimension-reduction plot style (UMAP, PCA, t-SNE).
#   source("assets/theme_embedding.R")
#   ggplot(df, aes(UMAP1, UMAP2)) + geom_point() + theme_embedding()
#
# Adds arrowed, corner-cropped axis guides and blanks the ticks/text — the
# convention for embeddings whose coordinates are unitless. Returns a list of
# ggplot2 layers, so it is *added* on top of the active theme (theme_nature());
# it does not replace it. Requires the `legendry` package for the fractional
# axis cap.
#
# Args:
#   axis  - axis guide; default crops the line to the first 20% of the panel
#           (legendry::guide_axis_base(cap = I(c(-Inf, 0.2)))).
#   arrow - grid::arrow() for the axis lines.
#
# Example corner label:
#   ... + theme_embedding() +
#     annotate("text", x = -Inf, y = -Inf, label = "UMAP",
#              hjust = -0.2, vjust = -0.6, size = 8 / .pt)
theme_embedding <- function(
    axis = legendry::guide_axis_base(cap = I(c(-Inf, 0.2))),
    arrow = grid::arrow(
      type = "closed",
      angle = 25,
      length = grid::unit(0.4, "lines")
    )) {
  list(
    ggplot2::theme(
      axis.ticks   = ggplot2::element_blank(),
      axis.text    = ggplot2::element_blank(),
      aspect.ratio = 1,
      axis.line    = ggplot2::element_line(arrow = arrow),
      # unitless orientation labels (PaCMAP 1/2) → quietest tier (fs_small), not full axis-title size
      axis.title   = ggplot2::element_text(hjust = 0, size = fs_small)
    ),
    ggplot2::guides(x = axis, y = axis)
  )
}
