library(tidyverse)
library(ggbrandon)
library(patchwork)
library(ggh4x)
library(legendry)
library(shadowtext)
library(ggrepel)

theme_set(theme_brandon(base_line_size = 0.3, base_size = 7))
# ----------------------------------------------------------------------
output_dir <- "figures/sf3"
source("make_figures/cfg.R")

branches <- c(
  "branch: 12,6,9,3,4,8",
  "branch: 12,6,9,3,10"
)
rep_factors <- c("Cebpa", "Fli1", "Bhlhe40", "Kdm1a")
rep_factors_pal <- setNames(
  RColorBrewer::brewer.pal(length(rep_factors), "Set2"),
  rep_factors
)
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

add_label <- function(df, x_col, y_col, priority, n_extra = 5) {
  x <- df[[x_col]]
  y <- df[[y_col]]
  x_n <- (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))
  y_n <- (y - min(y, na.rm = TRUE)) / diff(range(y, na.rm = TRUE))
  dm <- as.matrix(dist(cbind(x_n, y_n)))
  diag(dm) <- Inf
  nn_dist <- apply(dm, 1, min)
  is_priority <- df$factor %in% priority
  rank_order <- order(-nn_dist)
  show <- is_priority
  n_shown <- 0
  for (i in rank_order) {
    if (show[i]) next
    if (n_shown >= n_extra) break
    show[i] <- TRUE
    n_shown <- n_shown + 1
  }
  df %>% mutate(
    show_label = show,
    label_text = ifelse(show, factor, "")
  )
}

load_dat <- function(fname) {
  require(dplyr)
  require(tidyr)

  df_ps <-
    fname |>
    arrow::read_parquet() |>
    filter(
      branch %in% branches, # nolint
      is.na(subbranch) # nolint
    )

  dat_n <- df_ps |>
    filter(score < 0) |> # nolint
    group_by(branch, factor) |>
    summarise(score = sum(abs(score))) |> # nolint
    pivot_wider(names_from = branch, values_from = score) |>
    add_label("branch: 12,6,9,3,4,8", "branch: 12,6,9,3,10", rep_factors)

  dat_p <- df_ps |>
    filter(score > 0) |> # nolint
    group_by(branch, factor) |> # nolint
    summarise(score = sum(score)) |>
    pivot_wider(names_from = branch, values_from = score) |>
    add_label("branch: 12,6,9,3,4,8", "branch: 12,6,9,3,10", rep_factors)

  list(dat_n = dat_n, dat_p = dat_p)
}

dat <- load_dat("box/results/pgd_celloracle/perturb_scores.parquet")
dat2 <-
  load_dat("box/results/pgd_celloracle_overexpress/perturb_scores.parquet")

dat_n <- dat$dat_n
dat_p <- dat$dat_p
dat2_n <- dat2$dat_n
dat2_p <- dat2$dat_p

lim_n <- c(0, max(c(dat_n$`branch: 12,6,9,3,4,8`, dat_n$`branch: 12,6,9,3,10`), na.rm = TRUE))
lim_p <- c(0, max(c(dat_p$`branch: 12,6,9,3,4,8`, dat_p$`branch: 12,6,9,3,10`), na.rm = TRUE))
lim2_n <- c(0, max(c(dat2_n$`branch: 12,6,9,3,4,8`, dat2_n$`branch: 12,6,9,3,10`), na.rm = TRUE))
lim2_p <- c(0, max(c(dat2_p$`branch: 12,6,9,3,4,8`, dat2_p$`branch: 12,6,9,3,10`), na.rm = TRUE))

p1 <- dat_n |>
  ggplot(aes(
    x = `branch: 12,6,9,3,4,8`,
    y = `branch: 12,6,9,3,10`
  )) +
  geom_point(
    data = dat_n |> filter(!factor %in% rep_factors),
    fill = "grey50",
    shape = 21,
    color = "white",
    size = 1.2,
    stroke = 0.2
  ) +
  geom_point(
    aes(color = rep_factors_pal[factor]),
    dat_n |> filter(factor %in% rep_factors),
    shape = 16,
    size = 5,
    alpha = 0.3
  ) +
  geom_point(
    aes(fill = rep_factors_pal[factor]),
    dat_n |> filter(factor %in% rep_factors),
    shape = 21,
    color = "white",
    size = 1.2,
    stroke = 0.2
  ) +
  geom_text_repel(
    aes(
      label = label_text,
      fontface = ifelse(factor %in% rep_factors, "bold.italic", "italic")
    ),
    size = 5 / .pt,
    segment.size = 0.25,
    max.overlaps = Inf,
    min.segment.length = 0,
    max.iter = 1e5
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  labs(
    x = "branch 3: 12,6,9,3,4,8: negative PS sum",
    y = "branch 2: 12,6,9,3,10:\nnegative PS sum"
  ) +
  coord_fixed(clip = "off", xlim = lim_n, ylim = lim_n) +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(
      size = 6,
      color = RColorBrewer::brewer.pal(8, "PiYG")[1]
    )
  )

p2 <- dat_p |>
  ggplot(aes(
    x = `branch: 12,6,9,3,4,8`,
    y = `branch: 12,6,9,3,10`
  )) +
  geom_point(
    data = dat_p |> filter(!factor %in% rep_factors),
    fill = "grey50",
    shape = 21,
    color = "white",
    size = 1.2,
    stroke = 0.2
  ) +
  geom_point(
    aes(color = rep_factors_pal[factor]),
    dat_p |> filter(factor %in% rep_factors),
    shape = 16,
    size = 5,
    alpha = 0.3
  ) +
  geom_point(
    aes(fill = rep_factors_pal[factor]),
    dat_p |> filter(factor %in% rep_factors),
    shape = 21,
    color = "white",
    size = 1.2,
    stroke = 0.2
  ) +
  geom_text_repel(
    aes(
      label = label_text,
      fontface = ifelse(factor %in% rep_factors, "bold.italic", "italic")
    ),
    size = 5 / .pt,
    segment.size = 0.25,
    max.overlaps = Inf,
    min.segment.length = 0,
    max.iter = 1e5
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  labs(
    x = "branch 3: 12,6,9,3,4,8: positive PS sum",
    y = "branch 2: 12,6,9,3,10:\npositive PS sum"
  ) +
  coord_fixed(clip = "off", xlim = lim_p, ylim = lim_p) +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(
      size = 6,
      color = RColorBrewer::brewer.pal(8, "PiYG")[8]
    )
  )

p3 <- dat2_n |>
  ggplot(aes(
    x = `branch: 12,6,9,3,4,8`,
    y = `branch: 12,6,9,3,10`
  )) +
  geom_point(
    data = dat2_n |> filter(!factor %in% rep_factors),
    fill = "grey50",
    shape = 21,
    color = "white",
    size = 1.2,
    stroke = 0.2
  ) +
  geom_point(
    aes(color = rep_factors_pal[factor]),
    dat2_n |> filter(factor %in% rep_factors),
    shape = 16,
    size = 5,
    alpha = 0.3
  ) +
  geom_point(
    aes(fill = rep_factors_pal[factor]),
    dat2_n |> filter(factor %in% rep_factors),
    shape = 21,
    color = "white",
    size = 1.2,
    stroke = 0.2
  ) +
  geom_text_repel(
    aes(
      label = label_text,
      fontface = ifelse(factor %in% rep_factors, "bold.italic", "italic")
    ),
    size = 5 / .pt,
    segment.size = 0.25,
    max.overlaps = Inf,
    min.segment.length = 0,
    max.iter = 1e5
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  labs(
    x = "branch 3: 12,6,9,3,4,8: negative PS sum",
    y = "branch 2: 12,6,9,3,10:\nnegative PS sum"
  ) +
  coord_fixed(clip = "off", xlim = lim2_n, ylim = lim2_n) +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(
      size = 6,
      color = RColorBrewer::brewer.pal(8, "PiYG")[1]
    )
  )

p4 <- dat2_p |>
  ggplot(aes(
    x = `branch: 12,6,9,3,4,8`,
    y = `branch: 12,6,9,3,10`
  )) +
  geom_point(
    data = dat2_p |> filter(!factor %in% rep_factors),
    fill = "grey50",
    shape = 21,
    color = "white",
    size = 1.2,
    stroke = 0.2
  ) +
  geom_point(
    aes(color = rep_factors_pal[factor]),
    dat2_p |> filter(factor %in% rep_factors),
    shape = 16,
    size = 5,
    alpha = 0.3
  ) +
  geom_point(
    aes(fill = rep_factors_pal[factor]),
    dat2_p |> filter(factor %in% rep_factors),
    shape = 21,
    color = "white",
    size = 1.2,
    stroke = 0.2
  ) +
  geom_text_repel(
    aes(
      label = label_text,
      fontface = ifelse(factor %in% rep_factors, "bold.italic", "italic")
    ),
    nudge_y = ifelse(dat2_p$factor == "Fli1", -15, 0),
    size = 5 / .pt,
    segment.size = 0.25,
    max.overlaps = Inf,
    min.segment.length = 0,
    max.iter = 1e5
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  labs(
    x = "branch 3: 12,6,9,3,4,8: positive PS sum",
    y = "branch 2: 12,6,9,3,10:\npositive PS sum"
  ) +
  scale_x_continuous(breaks = c(0, 25, 50, 75)) +
  scale_y_continuous(breaks = c(0, 25, 50, 75)) +
  coord_fixed(clip = "off", xlim = lim2_p, ylim = lim2_p) +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(
      size = 6,
      color = RColorBrewer::brewer.pal(8, "PiYG")[8]
    )
  )


combined_ko <- p1 + plot_spacer() + p2 +
  plot_layout(nrow = 1, widths = c(1, 0.15, 1)) +
  plot_annotation(
    title = "KO analysis",
    theme = theme(
      plot.title = element_text(size = 7, margin = margin(t = 3), hjust = 0.5)
    )
  ) &
  theme(plot.margin = margin(0, 0, 0, 0))
save_figure(output_dir, "scatterplots_ps_ko",
  plot = combined_ko, width = 4.5, height = 2.3
)

combined_oe <- p3 + plot_spacer() + p4 +
  plot_layout(nrow = 1, widths = c(1, 0.15, 1)) +
  plot_annotation(
    title = "Overexpression analysis",
    theme = theme(
      plot.title = element_text(size = 7, margin = margin(t = 3), hjust = 0.5)
    )
  ) &
  theme(plot.margin = margin(0, 0, 0, 0))
save_figure(output_dir, "scatterplots_ps_overexpression",
  plot = combined_oe, width = 4.5, height = 2.3
)
