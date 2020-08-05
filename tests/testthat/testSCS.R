if(R.Version()$minor == "5.3"){
  test_that("Whether stim_cell_selector selects the same number of cells.",{
    selected_data <- stim_cell_selector(chi11_1k$expr_data, chi11_1k$state_markers,
                                        chi11_1k$cluster_col, chi11_1k$stim_label,
                                        chi11_1k$unstim_label, seed_val = 123,
                                        umap = TRUE, umap_cells = 50, verbose = TRUE)
    expect_equal(nrow(selected_data$selected_expr_data), 16055)
  })
} else if(R.Version()$minor == "6.3"){
  test_that("Whether stim_cell_selector selects the same number of cells.",{
    selected_data <- stim_cell_selector(chi11_1k$expr_data, chi11_1k$state_markers,
                                        chi11_1k$cluster_col, chi11_1k$stim_label,
                                        chi11_1k$unstim_label, seed_val = 123,
                                        umap = TRUE, umap_cells = 50, verbose = TRUE)
    #expect_equal(nrow(selected_data$selected_expr_data), 16023)
    expect_equal(nrow(selected_data$selected_expr_data), 17218)

  })
}

if(R.Version()$minor == "6.3"){
  test_that("Whether stim_cell_selector_single_marker selects the same number of cells.",{
    selected_data <- stim_cell_selector_single_marker(chi11_1k$expr_data, chi11_1k$state_markers,
                                                      chi11_1k$cluster_col, chi11_1k$stim_label,
                                                      chi11_1k$unstim_label, seed_val = 123,
                                                      verbose = TRUE)
    expect_equal(nrow(selected_data$selected_expr_data), 58727)
  })
}
