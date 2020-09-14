#' @title Stim Cell Selector Single Marker
#' @description Function to select cells from the stimulated samples that
#'              have likely responded to the stimulant.
#'
#' @param dat              A tibble with the single cell data. Cells on rows
#'                         and variables/markers on columns.
#' @param state_markers    A character vector with the labels of state
#'                         markers from the stimulation panel.
#' @param cluster_col      Column in the tibble with the cluster IDs.
#' @param stim_lab         A character of stim label(s).
#' @param unstim_lab       A character of unstim label(s).
#' @param seed_val         Seed value (integer) for \code{\link{kmeans}} clustering.
#'                         Default is NULL for no seed value.
#' @param verbose          Logical. To make function more verbose. Default is FALSE.
#'
#' @return A list with tibbles for expression data for the selected cells,
#'         data to plot stacked bar plots, data to plot UMAP plots, and
#'         parameters passed to the function.
#' @importFrom tibble as_tibble
#' @import dplyr ggplot2
#' @importFrom stats fisher.test kmeans median
#' @export
#'
#' @examples
#' selected_data <- stim_cell_selector_single_marker(chi11_1k$expr_data, chi11_1k$state_markers,
#'                   chi11_1k$cluster_col, chi11_1k$stim_label,
#'                   chi11_1k$unstim_label, seed_val = 123, verbose = FALSE)
stim_cell_selector_single_marker <- function(dat, state_markers, cluster_col, stim_lab, unstim_lab,
                               seed_val = NULL, verbose = FALSE){
  # For debugging.
  # verbose <- TRUE
  # dat <- chi11_1k$expr_data
  # state_markers <- chi11_1k$state_markers[1]
  # cluster_col <- chi11_1k$cluster_col
  # stim_lab <- chi11_1k$stim_label[1]
  # unstim_lab <- chi11_1k$unstim_label
  # seed_val <- 123

  # Check argument accuracy.
  if(verbose == TRUE){message(paste("Selection process started on", date()))}

  # Standardize column names and state marker labels.
  if(verbose == TRUE){message("Relacing - with _ in state marker labels and expression data column names.")}
  state_markers <- gsub("-", "_", state_markers)
  cols <- colnames(dat)
  cols <- gsub("-", "_", cols)
  colnames(dat) <- cols
  cols_not_state <- setdiff(cols, state_markers)

  # Fetch cluster labels from the expression data frame.
  clusters <- as.character(unique(dat[[cluster_col]]))

  total_combinations <- length(clusters) * length(stim_lab) * length(state_markers)
  if(verbose == TRUE){message(paste0("Total combinations to go through: ", length(clusters), " clusters X ", length(stim_lab),
                                     " stimulation types X ", length(state_markers), " state markers = ", total_combinations,"."))}

  # Initialize empty output data frames.
  df_stim_out <- data.frame(matrix(ncol = length(cols_not_state) + 2, nrow = 0))
  df_summary_out <- data.frame(matrix(ncol = 6, nrow = 0))
  stacked_bar_plot_data <- data.frame(matrix(ncol = 9, nrow = 0))

  # Set counter for the number of combinations processed.
  counter <- 0

  # Main nested for loop.
  for(stim in stim_lab){
    # stim <- stim_lab[1]
    for(cluster in clusters){
      # cluster <- clusters[1]
      for(state in state_markers){
        # state <- state_markers[4]

        # Set seed value for the k-means clustering.
        set.seed(seed_val)
        counter <- counter + 1
        if(verbose == TRUE){message(paste0("## Combination ", counter,"/",total_combinations,"."))}

        # Initialize an empty temporary output data frame.
        df_out_temp <- data.frame(matrix(ncol = length(cols), nrow = 0))

        # For the selected meta cluster, select cells from the selected stim type and un-stimulated cells.
        dat_stim_unstim <- dat[dat$stim_type == stim | dat$stim_type %in% unstim_lab, ]
        dat_stim_unstim <- dat_stim_unstim[dat_stim_unstim[,cluster_col] == cluster,]

        # Identify the index of the stimulated and unstimulated cells in the selected data frame.
        stim_idx <- dat_stim_unstim$stim_type == stim
        unstim_idx <- dat_stim_unstim$stim_type %in% unstim_lab

        # From dat_stim_unstim data frame select expression values for all the state markers
        # and carry out k-means (k = 2) clustering.
        if(verbose == TRUE){message(paste("Carrying out k-means clustering on cells from", cluster, "-", stim, "-", state, "+ unstim."))}
        dat_state <- dat_stim_unstim[, state]
        k_results <- kmeans(dat_state, 2)

        # Get cluster IDs for stim and unstim cells.
        clust_stim <- k_results$cluster[stim_idx]
        clust_unstim <- k_results$cluster[unstim_idx]

        # Fisher's Exact Test.
        if(verbose == TRUE){message("Carrying out Fisher's exact test.")}

        # Create a contingency table and carry out F-test.
        # Count the number of stimulated and un-stimulated cells that belong to cluster 1 and 2.
        total_stim_count <- as.numeric(nrow(dat_stim_unstim[stim_idx,]))
        total_unstim_count <- as.numeric(nrow(dat_stim_unstim[unstim_idx,]))

        stim_1 <- as.numeric(length(clust_stim[clust_stim == 1]))
        stim_1_ <- stim_1/total_stim_count
        stim_2 <- as.numeric(length(clust_stim[clust_stim == 2]))
        stim_2_ <- stim_2/total_stim_count

        unstim_1 <- as.numeric(length(clust_unstim[clust_unstim == 1]))
        unstim_1_ <- unstim_1/total_unstim_count
        unstim_2 <- as.numeric(length(clust_unstim[clust_unstim == 2]))
        unstim_2_ <- unstim_2/total_unstim_count

        con_tab <-  matrix(c(round(stim_1_ * 100), round(stim_2_ * 100), round(unstim_1_ * 100), round(unstim_2_ * 100)), nrow = 2, ncol = 2, dimnames = list(c("Cluster1", "Cluster2"), c("Stim", "Unstim")))

        message(con_tab)
        message(unstim_lab)
        f_test <- fisher.test(con_tab) # Fisher's exact test.
        f_p_val <- f_test$p.value

        # Select cells and save data only for the combinations that pass
        # F-test.
        if(f_p_val < 0.05){
          if(verbose == TRUE){message("Fisher's exact test significant.")}

          # Identify responding stimulated cells cluster.
          clust_1_ratio <- stim_1_ / unstim_1_
          clust_2_ratio <- stim_2_ / unstim_2_
          if(clust_1_ratio > clust_2_ratio){
            unstim_clust = 2
            stim_clust = 1
          } else {
            unstim_clust = 1
            stim_clust = 2
          }

          # Fetch responding stimulated cells from the indentified cluster.
          stim_cells_resp <- dat_stim_unstim[k_results$cluster == stim_clust,]
          stim_cells_resp <- stim_cells_resp[stim_cells_resp$stim_type == stim, ]

          # Select median marker expression for each state marker for the selected stim cells.
          tm <- stim_cells_resp[stim_cells_resp$stim_type == stim, state]
          med_exp <- median(tm[[1]])


          # tm <- dat_stim_unstim[k_results$cluster == stim_clust,]
          # tm <- tm[tm$stim_type == stim, state]
          # med_exp <- median(tm[[1]])

          # Calculate fold change from non-responsive to responsive state for stim cells.
          if(stim_clust == 1){
            fc <- round((stim_1/stim_2), digits = 2)
          }else if(stim_clust == 2){
            fc <- round((stim_2/stim_1), digits = 2)
          }
          # if(stim_clust == 1){
          #   fc <- round((stim_1 - stim_2)/stim_1, digits = 2)
          # }else if(stim_clust == 2){
          #   fc <- round((stim_2 - stim_1)/stim_2, digits = 2)
          # }

          # Selected stimulated cells.
          # stim_cells <- dat_stim_unstim[k_results$cluster == stim_clust,]
          # stim_cells <- stim_cells[as.character(stim_cells$stim_type) == stim, ]

          if(dim(stim_cells_resp)[1] != 0){
            # Combine data in a data frame and also generate a summary table.
            stim_cells_temp <- cbind(stim_cells_resp[cols_not_state], "state_marker" =  state, stim_cells_resp[state])
            names(stim_cells_temp)[names(stim_cells_temp) == state] <- "state_marker_exp"
            df_stim_out <- rbind(df_stim_out, stim_cells_temp)

            summary_table <- cbind(cluster = cluster, stim_type = stim, state_marker = state,
                                   f_p_value = f_p_val, fold_change = fc, state_maker_med_exp = med_exp)
            df_summary_out <- rbind(df_summary_out, summary_table)

            # Stacked + percent bar plots with stim and unstim on the x-axis and cell counts
            # in k-means clusters on the y-axis.
            stacked_data <- data.frame(stim_status = c("unstim", "unstim", "stim", "stim"),
                                       k_cluster = c("cluster1", "cluster2", "cluster1", "cluster2"),
                                       count = c(unstim_1_, unstim_2_, stim_1_, stim_2_))

            stacked_temp <- cbind("cluster" = cluster, "stim_type" = stim, "state_marker" = state, "f_p_val" = f_p_val,
                                "stim_clust" = stim_clust, "fold_change" = fc, stacked_data)
            stacked_bar_plot_data <- rbind(stacked_bar_plot_data, stacked_temp)
          }
        } else {
          if(verbose == TRUE){message("Fisher' exact test non-sgnificant. Skipping further steps.")}
        }
      }
    }
  }

  return_list <- list("selected_expr_data" = as_tibble(df_stim_out), "summary" = as_tibble(df_summary_out),
                      "stacked_bar_plot_data" = as_tibble(stacked_bar_plot_data),
                      "state_markers" = state_markers, "cluster_col" = cluster_col,
                      "stim_label" = stim_lab, "unstim_label" = unstim_lab,
                      "seed_val" = seed_val)

  if(verbose == TRUE){message(paste("Selection process finished on", date()))}
  return(return_list)
}


