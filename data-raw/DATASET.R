## code to prepare `DATASET` dataset goes here
library(tidyverse)

# Prepare CyTOF dataset.
panel <- read_tsv(file.path("data-raw", "phospho_panel_v3.txt"), show_col_types = FALSE)
type_markers <- panel[panel$marker_class == "type", ][["antigen"]]
type_markers <- gsub("-", "_", type_markers)
state_markers <- panel[panel$marker_class == "state", ][["antigen"]]
state_markers <- gsub("-", "_", state_markers)

stim_label <- c("A", "T", "L", "G")
unstim_label <- "U"

cluster_col <- "merging1"

chi11_1k_expr <- read_tsv(file.path("data-raw", "chi11_1k.tsv"), show_col_types = FALSE)
set.seed(12345)
chi11_1k_expr <- group_by(chi11_1k_expr, stim_type, merging1) %>%
  mutate("ncount" = n()) %>%
  dplyr::filter(ncount > 500) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  droplevels()

chi11 <- list("expr_data" = chi11_1k_expr, "type_markers" = type_markers, "state_markers" = state_markers,
                 "cluster_col" = cluster_col, "stim_label" = stim_label, "unstim_label" = unstim_label)

usethis::use_data(chi11, compress = "xz", overwrite = TRUE)
