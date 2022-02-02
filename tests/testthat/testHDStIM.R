if(R.Version()$major == "3" & R.Version()$minor == "6.3"){
  test_that("Whether HDStIM maps the same number of cells.",{
    mapped_data <- HDStIM(chi11_1k$expr_data, chi11_1k$state_markers,
                                        chi11_1k$cluster_col, chi11_1k$stim_label,
                                        chi11_1k$unstim_label, seed_val = 123,
                                        umap = TRUE, umap_cells = 50,
                                        verbose = FALSE)
    #expect_equal(nrow(mapped_data$selected_expr_data), 16023)
    resp_data <- filter(mapped_data$response_mapping_main, stim_type != "U" & k_cluster_id == responding_cluster)
    expect_equal(nrow(resp_data), 17204)

  })
}

if(R.Version()$major == "4"){
  test_that("Whether HDStIM maps the same number of cells.",{
    mapped_data <- HDStIM(chi11_1k$expr_data, chi11_1k$state_markers,
                                        chi11_1k$cluster_col, chi11_1k$stim_label,
                                        chi11_1k$unstim_label, seed_val = 123,
                                        umap = TRUE, umap_cells = 50,
                                        verbose = FALSE)
    #expect_equal(nrow(mapped_data$selected_expr_data), 16023)
    resp_data <- filter(mapped_data$response_mapping_main, stim_type != "U" & k_cluster_id == responding_cluster)
    expect_equal(nrow(resp_data), 17204)

  })
}
