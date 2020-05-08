#' @title Stim Cell Selector
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
#' @param umap             Boolean (T/F) to carry out UMAP on the selected cells.
#'                         Default is FALSE to skip UMAP calculation.
#' @param umap_cells       An integer; for calculating UMAPs take a minimum of \code{umap_cells} per
#'                         cluster or the total number of cells if the cluster size
#'                         is smaller than \code{umap_cells}. Default is NULL.
#'
#' @return A list with tibbles for expression data for the selected cells,
#'         data to plot stacked bar plots, and data to plot UMAP plots.
#' @importFrom tibble as_tibble
#' @import dplyr ggplot2
#' @importFrom stats fisher.test kmeans median
#' @importFrom uwot umap
#' @export
#'
#' @examples
#' selected_data <- stim_cell_selector(chi11_1k$expr_data, chi11_1k$state_markers,
#'                   chi11_1k$cluster_col, chi11_1k$stim_label,
#'                   chi11_1k$unstim_label, seed_val = 123, umap = FALSE, umap_cells = NULL)
#'
#' @note To reduce verbosity use \code{suppressMessages(stim_cells_selector())}.
stim_cell_selector <- function(dat, state_markers, cluster_col, stim_lab, unstim_lab, seed_val = NULL, umap = F, umap_cells = NULL){
  # For debugging.
  # dat <- chi11_1k$expr_data
  # state_markers <- chi11_1k$state_markers
  # cluster_col <- chi11_1k$cluster_col
  # stim_lab <- chi11_1k$stim_label
  # unstim_lab <- chi11_1k$unstim_label
  # seed_val <- 123
  # umap <- TRUE
  # umap_cells <- 50

  # Check argument accuracy.
  if(umap == TRUE & is.null(umap_cells)){
    stop("If umap is set to TRUE then please pass a value for umap_cells.")
  }

  message(paste("Selection process started on", date()))

  # Standardize column names and state marker labels.
  state_markers <- gsub("-", "_", state_markers)
  cols <- colnames(dat)
  cols <- gsub("-", "_", cols)
  colnames(dat) <- cols

  # Fetch cluster labels from the expression data frame.
  clusters <- as.character(unique(dat[[cluster_col]]))

  total_combinations <- length(clusters) * length(stim_lab)
  message(paste0("Total combinations to go through: ", length(clusters), " clusters X ", length(stim_lab),
                 " stimulation types = ", total_combinations,"."))

  # Initialize empty output data frames.
  df_out <- data.frame(matrix(ncol = length(cols), nrow = 0))
  df_summary_out <- as_tibble(data.frame(matrix(ncol = 3 + length(state_markers), nrow = 0)))
  stacked_bar_plot_data <- data.frame(matrix(ncol = 8, nrow = 0))
  if(umap){
    umap_plot_data <- data.frame(matrix(ncol = 7, nrow = 0))
  }

  # Set counter for the number of combinations processed.
  counter <- 0

  # Main nest for loop.
  for(stim in stim_lab){
    for(cluster in clusters){
      # Set seed value for the k-means clustering.
      set.seed(seed_val)
      counter <- counter + 1
      message(paste0("## Combination ", counter,"/",total_combinations,"."))
      # Initialize an empty temporary output data frame.
      df_out_temp <- data.frame(matrix(ncol = length(cols), nrow = 0))

      # For the selected meta cluster, select cells from the selected stim type and un-stimulated cells.
      dat_stim_unstim <- dat[dat$stim_type == stim | dat$stim_type == unstim_lab, ]
      dat_stim_unstim <- dat_stim_unstim[dat_stim_unstim[,cluster_col] == cluster,]

      # Identify the index of the stimulated and unstimulated cells in the selected data frame.
      stim_idx <- dat_stim_unstim$stim_type == stim
      unstim_idx <- dat_stim_unstim$stim_type == unstim_lab

      # From dat_stim_unstim data frame select expression values for all the state markers
      # and carry out k-means (k = 2) clustering.
      message(paste("Carrying out k-means clustering on cells from", cluster, "-", stim, "+ unstim."))
      dat_state <- dat_stim_unstim[, state_markers]
      k_results <- kmeans(dat_state, 2)

      # Get cluster IDs for stim and unstim cells.
      clust_stim <- k_results$cluster[stim_idx]
      clust_unstim <- k_results$cluster[unstim_idx]

      # Fisher's Exact Test.
      message("Carrying out Fisher's exact test.")

      # Create a contingency table and carry out F-test.
      # Count the number of stimulated and un-stimulated cells that belong to cluster 1 and 2.
      stim_1 <- as.numeric(length(clust_stim[clust_stim == 1]))
      stim_2 <- as.numeric(length(clust_stim[clust_stim == 2]))

      unstim_1 <- as.numeric(length(clust_unstim[clust_unstim == 1]))
      unstim_2 <- as.numeric(length(clust_unstim[clust_unstim == 2]))

      total_stim_count <- as.numeric(nrow(dat_stim_unstim[stim_idx,]))
      total_unstim_count <- as.numeric(nrow(dat_stim_unstim[unstim_idx,]))

      con_tab <-  matrix(c(round((stim_1/total_stim_count) * 100), round((stim_2/total_stim_count) * 100), round((unstim_1/total_unstim_count) * 100), round((unstim_2/total_unstim_count) * 100)), nrow = 2, ncol = 2, dimnames = list(c("Cluster1", "Cluster2"), c("Stim", "Unstim")))

      f_test <- fisher.test(con_tab) # Fisher's exact test.
      f_p_val <- f_test$p.value

      # Select cells and save data only for the combinations that pass
      # F-test.
      if(f_p_val < 0.05){
        message("Fisher's exact test significant.")
        # Identify un-stimulated cells cluster on the basis of the cluster that has a higher number of cells.
        if(unstim_1 > unstim_2){
          unstim_clust = 1
          stim_clust = 2
        } else {
          unstim_clust = 2
          stim_clust = 1
        }

        # Select median marker expression for each state marker for the selected stim cells.
        tm <- dat_stim_unstim[k_results$cluster == stim_clust,]
        tm <- tm[tm$stim_type == stim, state_markers]
        med_exp <- apply(tm, 2, median)
        med_exp <-  as.data.frame(t(med_exp))

        # Fold change from non-responsive to responsive state for stim cells.
        if(stim_clust == 1){
          fc <- round((stim_1 - stim_2)/stim_1, digits = 2)
        }else if(stim_clust == 2){
          fc <- round((stim_2 - stim_1)/stim_2, digits = 2)
        }

        # Selected stimulated cells.
        stim_cells <- dat_stim_unstim[k_results$cluster == stim_clust,]
        stim_cells <- stim_cells[as.character(stim_cells$stim_type) == stim, ]

        # Select stimulated cells that are clustered in un-stimulated cells cluster.
        # These are the cells from the simulated samples that likely didn't respond
        # to the stimulant.
        stim_cells_no_resp <- dat_stim_unstim[k_results$cluster == unstim_clust,]
        stim_cells_no_resp <- stim_cells_no_resp[as.character(stim_cells_no_resp$stim_type) == stim, ]

        # Also select unstimulated cells to combine them later with the
        # simulated cells for the UMAP.
        unstim_cells <- dat_stim_unstim[k_results$cluster == unstim_clust,]
        unstim_cells <- unstim_cells[as.character(unstim_cells$stim_type) == unstim_lab, ]

        # Combine data in a data frame and also generate a summary table.
        df_out <- rbind(df_out, stim_cells)

        summary_table <- cbind(cluster = cluster, stim_type = stim, f_p_value = f_p_val, fold_change = fc, med_exp)
        df_summary_out <- rbind(df_summary_out, summary_table)


        # Stacked + percent bar plots with stim and unstim on the x-axis and cell counts
        # in k-means clusters on the y-axis.
        stacked_data <- data.frame(stim_status = c("unstim", "unstim", "stim", "stim"),
                                   k_cluster = c("cluster1", "cluster2", "cluster1", "cluster2"),
                                   count = c(unstim_1, unstim_2, stim_1, stim_2))

        stacked_temp <- cbind("cluster" = cluster, "stim_type" = stim, "f_p_val" = f_p_val,
                              "stim_clust" = stim_clust, "fold_change" = fc, stacked_data)
        stacked_bar_plot_data <- rbind(stacked_bar_plot_data, stacked_temp)

        # Generate UMAP on combined stim and unstim cells and colour them
        # accordingly.
        if(umap == TRUE){
          stim_cells_label <- cbind(stim_cells, "type" = "stim_cells_resp")
          stim_cells_no_resp_label <- cbind(stim_cells_no_resp, "type" = "stim_cells_no_resp")
          unstim_cells_label <- cbind(unstim_cells, "type" = "unstim_cells")
          comb_stim_unstim <- as_tibble(rbind(stim_cells_label, stim_cells_no_resp_label, unstim_cells_label))
          tot_of_cells <- nrow(comb_stim_unstim)

          # Take a minimum of 100k cells per cluster or the total no. of cells if the cluster size is smaller than 100k.
          no_of_cells <- min(umap_cells, tot_of_cells)
          comb_stim_unstim <-  sample_n(comb_stim_unstim, no_of_cells, replace = FALSE)
          input_umap <- comb_stim_unstim[,state_markers]

          # Run UMAP.
          message(paste0("Running UMAP for ", cluster, " - ", stim,"."))
          u_res <- uwot::umap(input_umap)

          umap_temp <- cbind("cluster" = cluster, "stim_type" = stim,
                             "tot_of_cells" = tot_of_cells, "no_of_cells" = no_of_cells,
                             "UMAP1" = u_res[,1], "UMAP2" = u_res[,2],
                             "cell_type" = as.character(comb_stim_unstim$type))
          umap_plot_data <- rbind(umap_plot_data, umap_temp)
          umap_plot_data$UMAP1 <- as.numeric(umap_plot_data$UMAP1)
          umap_plot_data$UMAP2 <- as.numeric(umap_plot_data$UMAP2)
        }

      } else {
        message("Fisher' exact test non-sgnificant. Skipping further steps.")
      }
    }
  }

  if(umap){
    return_list <- list("selected_expr_data" = as_tibble(df_out), "summary" = as_tibble(df_summary_out),
                        "stacked_bar_plot_data" = as_tibble(stacked_bar_plot_data),
                        "umap_plot_data" = as_tibble(umap_plot_data))
  }else{
    return_list <- list("selected_expr_data" = as_tibble(df_out), "summary" = as_tibble(df_summary_out),
                        "stacked_bar_plot_data" = as_tibble(stacked_bar_plot_data))
  }
 message(paste("Selection process finished on", date()))
 return(return_list)
}


