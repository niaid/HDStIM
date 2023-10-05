if(R.Version()$major == "4"){
  test_that("Test Boruta output.",{
    mapped_data <- HDStIM(chi11$expr_data, chi11$state_markers,
                                        chi11$cluster_col, chi11$stim_label,
                                        chi11$unstim_label, seed_val = 123,
                                        umap = FALSE)

    marker_ranking <- marker_ranking_boruta(mapped_data, path = NULL,
                                             n_cells = 1000, max_runs = 20, seed_val = 123,
                                             verbose = FALSE)

    pht <- plot_marker_ranking_heatmap(marker_ranking)

    expect_type(marker_ranking, "list")
    expect_type(pht, "S4")

  })
}
