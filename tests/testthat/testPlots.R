test_that("If plot_sbp returns 0.",{
  mapped_data <- HDStIM(chi11_1k$expr_data, chi11_1k$state_markers,
                                      chi11_1k$cluster_col, chi11_1k$stim_label,
                                      chi11_1k$unstim_label, seed_val = 123,
                                      umap = FALSE)

  s <- plot_K_Fisher(mapped_data, path = NULL, verbose = TRUE)
  expect_equal(s, 0)
})

test_that("If plot_umap returns 0.",{
  mapped_data <- HDStIM(chi11_1k$expr_data, chi11_1k$state_markers,
                                      chi11_1k$cluster_col, chi11_1k$stim_label,
                                      chi11_1k$unstim_label, seed_val = 123,
                                      umap = TRUE, umap_cells = 50)

  u <- plot_umap(mapped_data, path = NULL, verbose = TRUE)
  expect_equal(u, 0)
})

test_that("If plot_kde returns 0.",{
  mapped_data <- HDStIM(chi11_1k$expr_data, chi11_1k$state_markers,
                                      chi11_1k$cluster_col, chi11_1k$stim_label,
                                      chi11_1k$unstim_label, seed_val = 123,
                                      umap = FALSE)

  k <- plot_exprs(mapped_data, path = NULL, verbose = TRUE)
  expect_equal(k, 0)
})
