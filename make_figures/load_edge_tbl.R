library(tidyverse)


load_edge_tbl <- function(tree, df) {
  require(tibble)
  require(dplyr)
  require(stringr)
  require(purrr)
  require(tidyr)

  tree <- setNames(
    tree$detection.rate$detection.rate,
    rownames(tree$detection.rate)
  )

  # Helper to parse names into numeric vectors
  parse_nodes <- function(x) {
    if (grepl("^c\\(", x)) {
      str_extract_all(x, "\\d+")[[1]] %>% as.integer()
    } else {
      str_split(x, ":")[[1]] %>% as.integer()
    }
  }

  # Flattened version
  edge_tbl <- enframe(tree, name = "nodes", value = "detection.rate") %>%
    mutate(node_list = map(nodes, parse_nodes)) %>% # nolint
    select(-nodes) %>%
    mutate(edges = map2(node_list, detection.rate, ~ tibble( # nolint
      source = .x[-length(.x)],
      target = .x[-1],
      detection.rate = .y
    ))) %>%
    select(edges) %>% # nolint
    unnest(edges)

  map_umap_1 <- df %>%
    group_by(clusterid) %>% # nolint
    summarise_if(is.numeric, median) %>%
    select(clusterid, umap_1) %>% # nolint
    deframe()

  map_umap_2 <- df %>%
    group_by(clusterid) %>% # nolint
    summarise_if(is.numeric, median) %>%
    select(clusterid, umap_2) %>% # nolint
    deframe()

  edge_tbl %>%
    mutate(
      source = factor(source, levels(df$clusterid)),
      target = factor(target, levels(df$clusterid)),
      x = map_umap_1[source],
      y = map_umap_2[source],
      xend = map_umap_1[target], # nolint
      yend = map_umap_2[target]
    )
}
