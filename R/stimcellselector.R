#' @title Stim Cell Selector
#' @description Function to select cells from the stimulated samples that
#'              have likely responded to the stimulant.
#'
#' @param dat              A tibble with the single cell data. Cells on rows
#'                         and variables/markers on columns.
#' @param state_markers    A character vector with the labels of state
#'                         markers from the stimulation panel.
#' @param cluster_col      Column in the tibble with the cluster IDs.
#' @param stim_lab         A character vector of stim label(s).
#' @param unstim_lab       A character of unstim label(s).
#' @param seed_val         Seed value (integer) for \code{\link{kmeans}} clustering.
#'                         Default is NULL for no seed value.
#' @param umap             Boolean (T/F) to carry out UMAP on the selected cells.
#'                         Default is FALSE to skip UMAP calculation.
#' @param umap_cells       An integer; for calculating UMAPs take a minimum of \code{umap_cells} per
#'                         cluster or the total number of cells if the cluster size
#'                         is smaller than \code{umap_cells}. Default is NULL.
#' @param lr               Boolean (T/F) to carry out logistic regression.
#' @param lr_max_it        Maximum iteration (default = 50) for logistic regression.
#' @param verbose          Logical. To make function more verbose. Default is FALSE.
#'
#' @return A list with tibbles for expression data for the selected cells,
#'         data to plot stacked bar plots, data to plot UMAP plots, and
#'         parameters passed to the function.
#' @importFrom tibble as_tibble
#' @import dplyr ggplot2
#' @importFrom stats fisher.test kmeans median glm coefficients
#' @import skmeans
#' @importFrom uwot umap
#' @export
#'
#' @examples
#' selected_data <- stim_cell_selector(chi11_1k$expr_data, chi11_1k$state_markers,
#'                   chi11_1k$cluster_col, chi11_1k$stim_label,
#'                   chi11_1k$unstim_label, seed_val = 123, umap = FALSE, umap_cells = NULL,
#'                   lr = FALSE, lr_max_it = 50, verbose = FALSE)
#'
stim_cell_selector <- function(dat, state_markers, cluster_col, stim_lab, unstim_lab,
                               seed_val = NULL, umap = FALSE, umap_cells = NULL,
                               lr = FALSE, lr_max_it = 50, verbose = FALSE){
  # For debugging.
  # library(stimcellselector)
  # library(tidyverse)
  # library(Boruta)
  # dat <- chi11_1k$expr_data
  # state_markers <- chi11_1k$state_markers
  # cluster_col <- chi11_1k$cluster_col
  # stim_lab <- chi11_1k$stim_label
  # unstim_lab <- chi11_1k$unstim_label
  # seed_val <- 123
  # umap <- TRUE
  # umap_cells <- 50
  # verbose <- TRUE
  # lr <- FALSE
  # lr_max_it <- 50

  # dat <- dat[which(rowSums(dat[state_markers]) != 0),]

  # Check argument accuracy.
  if(umap == TRUE & is.null(umap_cells)){
    stop("If umap is set to TRUE then please pass a value for umap_cells.")
  }

  if(verbose == TRUE){message(paste("Selection process started on", date()))}

  # Standardize column names and state marker labels.
  if(verbose == TRUE){message("Relacing - with _ in state marker labels and expression data column names.")}
  state_markers <- gsub("-", "_", state_markers)
  cols <- colnames(dat)
  cols <- gsub("-", "_", cols)
  colnames(dat) <- cols

  # Fetch cluster labels from the expression data frame.
  clusters <- as.character(unique(dat[[cluster_col]]))

  total_combinations <- length(clusters) * length(stim_lab)
  if(verbose == TRUE){message(paste0("Total combinations to go through: ", length(clusters), " clusters X ", length(stim_lab),
                 " stimulation types = ", total_combinations,"."))}

  # Initialize empty output data frames.
  df_stim_out <- data.frame(matrix(ncol = length(cols), nrow = 0))
  df_summary_out <- as_tibble(data.frame(matrix(ncol = 3 + length(state_markers), nrow = 0)))
  stacked_bar_plot_data <- data.frame(matrix(ncol = 8, nrow = 0))
  if(umap){
    umap_plot_data <- data.frame(matrix(ncol = 7, nrow = 0))
  }
  df_f_fail_out <- data.frame(matrix(ncol = 2, nrow = 0))
  df_lr_out <- data.frame(matrix(nrow = 0, ncol = 8))
  df_all_f_out <- data.frame(matrix(nrow = 0, ncol = 7)) # May not be needed in the future.
  k_clust_data <- data.frame(matrix(nrow = 0 , ncol = 4 + length(state_markers)))

  # Set counter for the number of combinations processed.
  counter <- 0
  f_comb_no <- 0

  # Main nested for loop.
  for(stim in stim_lab){
    # stim <- stim_lab[1] # For debugging.
    for(cluster in clusters){
      # cluster <- clusters[2] # For debugging.
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
      if(verbose == TRUE){message(paste("Carrying out k-means clustering NF on cells from", cluster, "-", stim, "+ unstim."))}
      dat_state <- dat_stim_unstim[, state_markers]
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
      f_test <- fisher.test(con_tab) # Fisher's exact test.
      f_p_val <- f_test$p.value

      # Note p-value for all the Fisher's exact tests.
      df_all_f_out <- rbind(df_all_f_out,
                            data.frame("stim_type" = stim, "cluster" = cluster,
                                       "stim_clust1" = con_tab[1,1], "stim_clust2" = con_tab[2,1],
                                       "unstim_clust1" = con_tab[1,2], "unstim_clust2" = con_tab[2,2],
                                       "f_p_val" = f_p_val))

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

        # Combine responding stimulated cells in the output data frame.
        df_stim_out <- rbind(df_stim_out, stim_cells_resp)

        # Generate data for the summary table.
        # Select median marker expression for each state marker for the selected stim cells.
        tm <- stim_cells_resp[stim_cells_resp$stim_type == stim, state_markers]
        med_exp <- apply(tm, 2, median)
        med_exp <-  as.data.frame(t(med_exp))

        # Calculate fold change from non-responsive to responsive state for stim cells.
        if(stim_clust == 1){
          fc <- round((stim_1/stim_2), digits = 2)
        }else if(stim_clust == 2){
          fc <- round((stim_2/stim_1), digits = 2)
        }

        # Generate a summary table.
        summary_table <- cbind(cluster = cluster, stim_type = stim, f_p_value = f_p_val, fold_change = fc, med_exp)
        df_summary_out <- rbind(df_summary_out, summary_table)

        # Stacked + percent bar plots with stim and unstim on the x-axis and cell counts
        # in k-means clusters on the y-axis.
        stacked_data <- data.frame(stim_status = c("unstim", "unstim", "stim", "stim"),
                                   k_cluster = c("cluster1", "cluster2", "cluster1", "cluster2"),
                                   count = c(unstim_1_, unstim_2_, stim_1_, stim_2_))

        stacked_temp <- cbind("cluster" = cluster, "stim_type" = stim, "f_p_val" = f_p_val,
                              "stim_clust" = stim_clust, "fold_change" = fc, stacked_data)
        stacked_bar_plot_data <- rbind(stacked_bar_plot_data, stacked_temp)

        # Generate UMAP on combined stim and unstim cells and colour them
        # accordingly.
        if(umap == TRUE){
          # Select stimulated cells that are clustered in un-stimulated cells cluster.
          # These are the cells from the simulated samples that likely didn't respond
          # to the stimulant.
          stim_cells_no_resp <- dat_stim_unstim[k_results$cluster == unstim_clust,]
          stim_cells_no_resp <- stim_cells_no_resp[stim_cells_no_resp$stim_type == stim, ]

          # Fetch all the unstim cells.
          unstim_cells <- dat_stim_unstim[dat_stim_unstim$stim_type %in% unstim_lab, ]

          # Label cell types.
          stim_cells_label <- cbind(stim_cells_resp, "type" = "stim_cells_resp")
          stim_cells_no_resp_label <- cbind(stim_cells_no_resp, "type" = "stim_cells_no_resp")
          unstim_cells_label <- cbind(unstim_cells, "type" = "unstim_cells")
          comb_stim_unstim <- as_tibble(rbind(stim_cells_no_resp_label, stim_cells_label, unstim_cells_label))
          tot_of_cells <- nrow(comb_stim_unstim)

          # Take a minimum of "umap_cells" cells per cluster or the total no. of cells
          # if the cluster size is smaller than "umap_cells".
          no_of_cells <- min(umap_cells, tot_of_cells)
          comb_stim_unstim <-  sample_n(comb_stim_unstim, no_of_cells, replace = FALSE)
          input_umap <- comb_stim_unstim[,state_markers]

          # Run UMAP.
          if(verbose == TRUE){message(paste0("Running UMAP for ", cluster, " - ", stim,"."))}
          u_res <- uwot::umap(input_umap)

          umap_temp <- cbind("cluster" = cluster, "stim_type" = stim,
                             "tot_of_cells" = tot_of_cells, "no_of_cells" = no_of_cells,
                             "UMAP1" = u_res[,1], "UMAP2" = u_res[,2],
                             "cell_type" = as.character(comb_stim_unstim$type))
          umap_plot_data <- rbind(umap_plot_data, umap_temp)
          umap_plot_data$UMAP1 <- as.numeric(umap_plot_data$UMAP1)
          umap_plot_data$UMAP2 <- as.numeric(umap_plot_data$UMAP2)
        }

        # Logistic regression for each marker to test it's contribution in
        # k-means clustering.
        if(lr == TRUE){
          dat_for_lr <- mutate(dat_stim_unstim, k_cluster_id = k_results$cluster)
          dat_for_lr$k_cluster_id <- as.factor(dat_for_lr$k_cluster_id)
          for(marker in state_markers){
            # marker <- state_markers[1] # For Debugging.
            lr_form <- as.formula(paste0("k_cluster_id ~ ", marker))
            if(verbose == TRUE){message(paste0("Carring out logistic regression for ", marker))}
            lr_res <- suppressWarnings(glm(lr_form,family = "binomial", data = dat_for_lr, maxit = lr_max_it))
            est_marker <- coefficients(summary(lr_res))[2,1]
            st_er_marker <- coefficients(summary(lr_res))[2,2]
            z_marker <- coefficients(summary(lr_res))[2,3]
            pr_marker <- coefficients(summary(lr_res))[2,4]
            lr_aic <- lr_res$aic
            df_lr_out <- rbind(df_lr_out,
                               data.frame("stim_type" = stim, "cluster" = cluster, "state_marker" = marker,
                                          "lr_est" = est_marker, "lr_std_err" = st_er_marker,
                                          "lr_z_value" = z_marker, "lr_p_val" = pr_marker, "lr_aic" = lr_aic))
          }
        }

        # Generate data frame with statmarkers and k-means cluster identity per cell.
        f_comb_no <- f_comb_no + 1
        k_temp <- select(dat_stim_unstim, c("stim_type", all_of(cluster_col), all_of(state_markers))) %>%
          mutate("k_cluster_id" = k_results$cluster) %>%
          rename("cluster" = cluster_col) %>%
          mutate(k_cluster_id = as.factor(k_cluster_id)) %>%
          mutate("comb_no" = f_comb_no)
        k_clust_data <- rbind(k_clust_data, k_temp)

      } else {
        if(verbose == TRUE){message("Fisher' exact test non-sgnificant. Skipping further steps.")}
        df_f_fail_out <- rbind(df_f_fail_out, data.frame("stim_type" = stim, "cluster" = cluster))
      }
    }
  }

  # Generate return list.
  return_list <- list("selected_expr_data" = as_tibble(df_stim_out), "summary" = as_tibble(df_summary_out),
                      "stacked_bar_plot_data" = as_tibble(stacked_bar_plot_data),
                      "state_markers" = state_markers, "cluster_col" = cluster_col,
                      "stim_label" = stim_lab, "unstim_label" = unstim_lab,
                      "seed_val" = seed_val, "fisher_test_fail" = as_tibble(df_f_fail_out),
                      "all_fisher_p_val" = df_all_f_out, "k_clust_data" = k_clust_data)
  if(umap == TRUE){
    return_list[["umap_plot_data"]] <- as_tibble(umap_plot_data)
    return_list[["umap"]] <- umap
    return_list[["umap_cells"]] <- umap_cells
  }

  if(lr == TRUE){
    return_list[["lr_per_marker"]] <- as_tibble(df_lr_out)
  }

 if(verbose == TRUE){message(paste("Selection process finished on", date()))}
 return(return_list)
}


