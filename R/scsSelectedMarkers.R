#' @title SCS Selected Marker Pipeline
#' @description Meta function to run \code{stim_cell_selector} on
#'              the markers selected by \code{state_marker_selection_lmm}
#'              for each eligible stim type and cluster combination.
#'
#' @param dat                      A tibble with the single cell data. Cells on rows
#'                                 and variables/markers on columns.
#' @param selected_state_markers   Output data from \code{state_marker_selection_lmm} function.
#' @param cluster_col              Label for the column in the single cell data tibble that has the cluster information per cell.
#' @param stim_lab                 A character vector of stim label(s) for the \code{stim_cell_selector} function.
#' @param unstim_lab               A character of unstim label(s) for the \code{stim_cell_selector} function.
#' @param path                     Path to the directory where the output for this meta function will be saved.
#' @param anova_cutoff             FDR cutoff (Default is 0.05) to filter selected markers from \code{state_marker_selection_lmm} function.
#' @param seed_val                 Seed value (integer; default is NULL) for \code{\link{kmeans}} clustering in \code{stim_cell_selector_function}.
#' @param umap                     Boolean (T/F) to carry out UMAP on the selected cells. Default is FALSE to skip UMAP calculation.
#' @param umap_cells               An integer; for calculating UMAPs take a minimum of \code{umap_cells} per cluster
#'                                 or the total number of cells if the cluster size is smaller than \code{umap_cells}. Default is NULL.
#' @param verbose                  Logical. To make function more verbose. Default is FALSE.
#' @importFrom tibble as_tibble
#' @return A tibble with all the combinations for stim type
#'         and cluster that did not yield more than two
#'         state markers for the given anova fdr p-value cutoff.
#' @export
#' @examples
#' \donttest{
#' no_markers_comb <- scs_selected_markers(chi11_1k$expr_data, selected_state_markers,
#'                                         chi11_1k$cluster_col, chi11_1k$stim_label,
#'                                         chi11_1k$unstim_label, path = NULL,
#'                                         anova_cutoff = 0.05, seed_val = 123,
#'                                         umap = FALSE, umap_cells = NULL, verbose = FALSE)
#'}
scs_selected_markers <- function(dat, selected_state_markers, cluster_col, stim_lab, unstim_lab, path,
                                 anova_cutoff = 0.05, seed_val = NULL, umap = FALSE,
                                 umap_cells = NULL, verbose = FALSE){


  # For debugging.
  # dat <- chi11_1k$expr_data
  # anova_cutoff <- 0.05
  # cluster_col <- chi11_1k$cluster_col
  # stim_lab <- chi11_1k$stim_label
  # unstim_lab <- chi11_1k$unstim_label
  # seed_val <- 123
  # umap <- TRUE
  # umap_cells <- 50
  # verbose <- TRUE
  # path <- NULL

  # Check if path exists; if not then create it.
  if(!is.null(path)){
    if(!dir.exists(path)){
      if(verbose){message(paste("Creating %s folder", path))}
      dir.create(path, recursive = TRUE)
    } else {
      if(verbose){message(paste(path, "folder already exists. Output will be over written."))}
    }
  } else {message("WARNING: No path provided. Output will not be saved.")}

  # Get the selected marker data from the list.
  selected_markers_in <- selected_state_markers$df_fdr

  # Initiate an output data frame for combinations that did not yield any markers.
  df_no_markers <- data.frame(matrix(nrow = 0, ncol = 2))

  # Loop over stim and cluster combinatons and carry out SCS for each separately.
  clusters <- as.character(unique(dat[[cluster_col]]))
  for(stim in stim_lab){
    for(clust in clusters){
      selected_markers <- selected_markers_in[selected_markers_in$stim_type == stim & selected_markers_in$cluster == clust &
                                       selected_markers_in$anova_p_adj < anova_cutoff, ]
      # selected_markers <- selected_markers[selected_markers$fc > 1.2 | selected_markers$fc < 0.08,
      #                                      c("state_marker","anova_p_adj")]

      # Carry out SCS if there are atleast two selected markers
      # for the stim + clust combination.
      if(nrow(selected_markers) >= 2){
        selected_markers <- as.character(selected_markers$state_marker)
        if(verbose == TRUE){message(paste("No. of selected markers for stim", stim, "& cluster", clust, "=", length(selected_markers)))}

        dat_clust <- dat[dat[[cluster_col]] == clust, ]
        selected_data <-   stim_cell_selector(dat_clust, selected_markers, cluster_col, stim_lab = stim, unstim_lab,
                                              seed_val = seed_val, umap = umap, umap_cells = umap_cells,
                                              verbose = verbose)
        if(!is.null(path)){
          saveRDS(selected_data, file.path(path, paste0(stim, "_",clust, "_",anova_cutoff, ".rds")))
        }
      }else{
        df_no_markers <- rbind(df_no_markers, data.frame("stim_type" = stim, "cluster" = clust))
      }
    }
  }

  return(as_tibble(df_no_markers))
}
