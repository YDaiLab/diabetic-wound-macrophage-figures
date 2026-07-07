library(tidyverse)
source("assets/theme_nature.R")
library(patchwork)
library(ggh4x)
library(legendry)

theme_set(theme_nature())

# ----------------------------------------------------------------------
output_dir <- "figures/f2"
source("make_figures/cfg.R")
# ----------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- load_df()

df1 <- bind_rows(
  df %>%
    count(sample) %>%
    mutate(p = n / sum(n), c = 1),
  df %>%
    count(clusterid, sample) %>%
    group_by(clusterid) %>%
    mutate(p = n / sum(n), c = 2)
) %>%
  arrange(clusterid) %>%
  mutate(
    x = ifelse(is.na(clusterid), "Dataset", as.character(clusterid)),
    x = fct_inorder(x)
  )

df2 <- df1 %>%
  group_by(c, x) %>%
  summarise_at("n", sum) %>%
  mutate(
    nl = scales::comma(n),
    l = sprintf("%s\n(%s)", x, scales::comma(n))
  )

pa <- df2 %>%
  filter(c == 1) %>%
  ggplot(aes(x, n)) +
  geom_col(aes(fill = x)) +
  geom_text(
    aes(label = nl),
    vjust = -0.25,
    size = pt2mm(fs_small)
  ) +
  scale_fill_manual(values = cluster_pal) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 3),
    expand = expansion(mult = c(0, 0.2)),
    labels = scales::comma_format()
  ) +
  labs(
    x = NULL,
    y = "Count"
  ) +
  theme(
    axis.text.y = element_text(size = fs_small),
    strip.text = element_blank(),
    legend.position = "none"
  ) +
  coord_cartesian(clip = "off")

pb <- df2 %>%
  filter(c == 2) %>%
  ggplot(aes(x, n)) +
  geom_col(aes(fill = x)) +
  geom_text(
    aes(label = nl),
    vjust = -0.25,
    size = pt2mm(fs_small)
  ) +
  scale_fill_manual(values = cluster_pal) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 3),
    expand = expansion(mult = c(0, 0.2)),
    labels = scales::comma_format()
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = fs_small),
    strip.text = element_blank(),
    legend.position = "none"
  ) +
  coord_cartesian(clip = "off")

yintercept <- df %>%
  count(condition) %>%
  mutate(p = n / sum(n)) %>%
  filter(condition == "Diabetic") %>%
  pull(p)

pc <- df1 %>%
  filter(c == 1) %>%
  ggplot(aes(x, p)) +
  geom_col(aes(fill = sample)) +
  geom_hline(
    yintercept = yintercept,
    lty = 2,
    lwd = 0.3
  ) +
  scale_fill_manual(values = sample_pal) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    expand = expansion()
  ) +
  labs(
    x = NULL,
    y = "Proportion"
  ) +
  theme(
    axis.text.y = element_text(size = fs_small),
    strip.text = element_blank(),
    legend.position = "none"
  )

pd <- df1 %>%
  filter(c == 2) %>%
  ggplot(aes(x, p)) +
  geom_col(aes(fill = sample)) +
  geom_hline(
    yintercept = yintercept,
    lty = 2,
    lwd = 0.3
  ) +
  scale_fill_manual(values = sample_pal) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    expand = expansion()
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = fs_small),
    strip.text = element_blank(),
    legend.position = "none"
  )

design <- "
12
34
"
combined <- pa + pb + pc + pd +
  plot_layout(design = design, heights = c(1, 2), widths = c(1.5, 12))
save_figure(output_dir, "barplots", plot = combined, width = 3.5, height = 1.7)
