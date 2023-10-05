#' @title Marker ranking heatmap
#' @description
#' A consolidated heatmap showing the importance scores of all the state markers (X-axis) from
#' all the stimulation-cell population combinations that passed the Fisher's exact test (Y-axis).
#'
#' @param marker_ranking    Returned list from the \code{\link{marker_ranking_boruta}} function.
#' @import ComplexHeatmap circlize magrittr
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @return                 A ComplexHeatmap object
#' @export
#' @examples
#' \donttest{
#' mapped_data <- HDStIM(chi11$expr_data, chi11$state_markers,
#'                       chi11$cluster_col, chi11$stim_label,
#'                       chi11$unstim_label, seed_val = 123, umap = TRUE,
#'                       umap_cells = 50, verbose = FALSE)
#'
#' marker_ranking <- marker_ranking_boruta(mapped_data, path = NULL, n_cells = NULL,
#'                                         max_runs = 1000, seed_val = 123,
#'                                         verbose = 0)
#'
#' pht <- plot_marker_ranking_heatmap(marker_ranking)
#' }
plot_marker_ranking_heatmap <- function(marker_ranking){

  cell_population <- meanImp <- min_max <- state_marker <- stim_type <- NULL

  att_stats <- marker_ranking$attribute_stats %>%
    dplyr::select(stim_type, cell_population, state_marker, meanImp) %>%
    dplyr::group_by(stim_type, cell_population) %>%
    dplyr::mutate(min_max = (meanImp - min(meanImp)) /(max(meanImp)-min(meanImp))) %>%
    dplyr::select(stim_type, cell_population, state_marker, min_max) %>%
    dplyr::mutate("stim_pop" = paste0(stim_type, "::", cell_population)) %>%
    dplyr::ungroup()

  mat <- matrix(nrow = length(unique(att_stats$stim_pop)), ncol = length(unique(att_stats$state_marker)))
  colnames(mat) <- unique(att_stats$state_marker)
  rownames(mat) <- unique(att_stats$stim_pop)
  for(i in 1:nrow(att_stats)){
    rowi <- att_stats$stim_pop[i]
    coli <- as.character(att_stats$state_marker[i])
    val <- att_stats$min_max[i]
    mat[rowi, coli] <- val
  }

  sorted_cols <- order(colnames(mat))
  mat <- mat[, sorted_cols]

  sorted_rows <- order(rownames(mat))
  mat <- mat[sorted_rows,]

  r_breaks <- sapply(rownames(mat), function(x) strsplit(x, "::")[[1]][1]) %>%
    as.character()

  cleaned_row_names <- sapply(rownames(mat), function(x) strsplit(x, "::")[[1]][2]) %>%
    as.character()

  rownames(mat) <- cleaned_row_names

  col_fun = colorRamp2(c(0, 1), c("black", "yellow"))
  col_fun(seq(0, 1))

  hmap <- ComplexHeatmap::Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE,
                                  col = col_fun, row_split = r_breaks,
                                  heatmap_legend_param = list(title = "Importance"))

                                  # column_names_gp = gpar(fontsize = 10),
                                  # row_names_gp = gpar(fontsize = 10))

  # png(filename = file.path("/Users/farmerr2/sandbox/development", "marker_ranking_heatmap.png") , width = 7,
  #    height = 5, units = "in", pointsize = 12, bg = "white", res = 600)
  # draw(hmap)
  # dev.off()

  return(hmap)
}
