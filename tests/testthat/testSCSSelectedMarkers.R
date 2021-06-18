# test_that("If scs_selected_markers returns 43 combinations with not enough state markers for the high dimensional SCS.",{
#   selected_state_markers <- state_marker_selection_lmm(chi11_1k$expr_data, chi11_1k$state_markers,
#                                                        chi11_1k$cluster_col,chi11_1k$stim_label,
#                                                        chi11_1k$unstim_label, REML = FALSE, verbose = FALSE)
#
#   no_markers_comb <- scs_selected_markers(chi11_1k$expr_data, selected_state_markers, chi11_1k$cluster_col,
#                                           chi11_1k$stim_label, chi11_1k$unstim_label, path = NULL,
#                                           anova_cutoff = 0.05, seed_val = 123,
#                                           umap = FALSE, umap_cells = NULL,
#                                           verbose = FALSE)
#
#   expect_equal(nrow(no_markers_comb), 43)
# })
