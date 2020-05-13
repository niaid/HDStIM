test_that("If plot_sbp returns 0.",{
  selected_data <- stim_cell_selector(chi11_1k$expr_data, chi11_1k$state_markers,
                                      chi11_1k$cluster_col, chi11_1k$stim_label,
                                      chi11_1k$unstim_label, seed_val = 123,
                                      umap = FALSE)

  s <- plot_sbp(selected_data, path = NULL, verbose = TRUE)
  expect_equal(s, 0)
})

test_that("If plot_umap returns 0.",{
  selected_data <- stim_cell_selector(chi11_1k$expr_data, chi11_1k$state_markers,
                                      chi11_1k$cluster_col, chi11_1k$stim_label,
                                      chi11_1k$unstim_label, seed_val = 123,
                                      umap = TRUE, umap_cells = 50)

  u <- plot_umap(selected_data, path = NULL, verbose = TRUE)
  expect_equal(u, 0)
})

