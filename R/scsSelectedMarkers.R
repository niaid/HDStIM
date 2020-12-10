#' SCS Selected Marker Pipeline
#'
#' @param dat
#' @param selected_state_markers
#' @param cluster_col
#' @param stim_lab
#' @param unstim_lab
#' @param path
#' @param anova_cutoff
#' @param seed_val
#' @param umap
#' @param umap_cells
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
#' no_markers_comb <- scs_selected_markers(chi11_1k$expr_data, selected_state_markers, chi11_1k$cluster_col,
#'                                         chi11_1k$stim_label, chi11_1k$unstim_label, path = NULL,
#'                                         anova_cutoff = 0.05, seed_val = 123, umap = FALSE,
#'                                         umap_cells = NULL, verbose = FALSE)
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
  } else {warning("No path provided. Output will not be saved.")}

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

  return(df_no_markers)
}
