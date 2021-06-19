if(R.Version()$major == "4"){
  test_that("Test Boruta output.",{
    selected_data <- stim_cell_selector(chi11_1k$expr_data, chi11_1k$state_markers,
                                        chi11_1k$cluster_col, chi11_1k$stim_label,
                                        chi11_1k$unstim_label, seed_val = 123,
                                        umap = FALSE)

    attribute_stats <- marker_ranking_boruta(selected_data, path = NULL,
                                             max_runs = 20, seed_val = 123,
                                             verbose = FALSE)

    expect_equal(nrow(attribute_stats), 350)

  })
}
