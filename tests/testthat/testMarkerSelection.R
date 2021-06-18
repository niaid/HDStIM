# test_that("If the state marker selection returns the same length of tibble.",{
#   selected_state_markers <- state_marker_selection_lmm(chi11_1k$expr_data, chi11_1k$state_markers,
#                                                    chi11_1k$cluster_col,chi11_1k$stim_label,
#                                                    chi11_1k$unstim_label, REML = FALSE, verbose = FALSE)
#   expect_equal(nrow(selected_state_markers$df_fdr), 800)
# })
