# Nature-journal-style ggplot2 theme.
#   source("assets/theme_nature.R")
#   theme_set(theme_nature())
#
# Conventions: black (not grey) for all text and lines, no gridlines, sans-serif
# font, white background, ~0.5pt lines. Type follows a small, fixed TYPE SCALE
# (the Nature 5-7pt window + an 8pt Arial Black panel letter) rather than one uniform
# size: primary text (titles) one step above secondary text (tick labels), which
# is the conventional figure hierarchy and still "quiet" (only two sizes in use).

# ---- Type scale (pt) ------------------------------------------------------
# Semantic roles, NOT arbitrary sizes. The journal dictates the window, so these
# are the whole hierarchy. Panels reference fs_* directly for annotations; for
# geom_text()/geom_label() (whose `size` is in MM) wrap with pt2mm().
fs_label <- 8 # panel letter a,b,c — BOLD, upright; added in Inkscape, plot.tag matches
fs_base  <- 7 # PRIMARY: axis titles, legend titles, strip/facet titles, key on-panel text
fs_small <- 6 # SECONDARY: axis tick labels, legend text, most in-panel annotation
fs_min   <- 5 # FLOOR: dense/fine text only — never go below (production rejects < 5pt)

# Panel-letter face: a heavy sans that carries its own weight, so plot.tag uses
# family = ff_tag with face = "plain" (NOT face = "bold", which would synthesize
# a fake-bold on top of an already-black weight).
ff_tag   <- "Arial Black"

# geom_text()/geom_label() take `size` in millimetres; convert a point size first:
#   geom_text(size = pt2mm(fs_small))
pt2mm <- function(pt) pt / (72.27 / 25.4) # ggplot2's `.pt`; self-contained

theme_nature <- function(base_size = fs_base, base_family = "") {
  line_size <- 0.4 # pt; Nature target is ~0.5pt at final print size
  secondary <- base_size - 1 # tick labels / legend text one step below titles (fs_small when base=fs_base)
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text              = ggplot2::element_text(color = "black",
                                                size = base_size),
      plot.title        = ggplot2::element_text(color = "black",
                                                size = base_size,
                                                face = "plain", hjust = 0),
      plot.subtitle     = ggplot2::element_text(color = "black",
                                                size = secondary),
      plot.tag          = ggplot2::element_text(color = "black",
                                                family = ff_tag,
                                                size = fs_label,
                                                face = "plain", hjust = 0),
      axis.title        = ggplot2::element_text(color = "black",
                                                size = base_size),
      axis.text         = ggplot2::element_text(color = "black",
                                                size = secondary),
      legend.title      = ggplot2::element_text(color = "black",
                                                size = base_size),
      legend.text       = ggplot2::element_text(color = "black",
                                                size = secondary),
      strip.text        = ggplot2::element_text(color = "black",
                                                size = base_size),

      line              = ggplot2::element_line(color = "black",
                                                linewidth = line_size),
      axis.line         = ggplot2::element_line(color = "black",
                                                linewidth = line_size),
      axis.ticks        = ggplot2::element_line(color = "black",
                                                linewidth = line_size),

      panel.grid        = ggplot2::element_blank(),
      panel.background  = ggplot2::element_blank(),
      plot.background   = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.key        = ggplot2::element_blank(),
      strip.background  = ggplot2::element_blank()
    )
}
