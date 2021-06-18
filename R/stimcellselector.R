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
#' @param verbose          Logical. To make function more verbose. Default is FALSE.
#'
#' @return A list with tibbles for expression data for the selected cells,
#'         data to plot stacked bar plots, data to plot UMAP plots, and
#'         parameters passed to the function.
#' @importFrom tibble as_tibble
#' @importFrom tidyselect all_of
#' @importFrom broom tidy
#' @import dplyr ggplot2
#' @importFrom stats fisher.test kmeans median
#' @import skmeans
#' @importFrom uwot umap
#' @export
#'
#' @examples
#' selected_data <- stim_cell_selector(chi11_1k$expr_data, chi11_1k$state_markers,
#'                   chi11_1k$cluster_col, chi11_1k$stim_label,
#'                   chi11_1k$unstim_label, seed_val = 123, umap = FALSE, umap_cells = NULL,
#'                   verbose = FALSE)
#'
stim_cell_selector <- function(dat, state_markers, cluster_col, stim_lab, unstim_lab,
                               seed_val = NULL, umap = FALSE, umap_cells = NULL,
                               verbose = FALSE){
  # For debugging.
  # library(stimcellselector)
  # library(tidyverse)
  # dat <- chi11_1k$expr_data
  # state_markers <- chi11_1k$state_markers
  # cluster_col <- chi11_1k$cluster_col
  # stim_lab <- chi11_1k$stim_label
  # unstim_lab <- chi11_1k$unstim_label
  # seed_val <- 123
  # umap <- TRUE
  # umap_cells <- 50
  # verbose <- TRUE


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

  # Convert input data to a tibble.
  dat <- as_tibble(dat)

  # Fetch cluster labels from the expression data frame.
  clusters <- as.character(unique(dat[[cluster_col]]))

  total_combinations <- length(clusters) * length(stim_lab)
  if(verbose == TRUE){message(paste0("Total combinations to go through: ", length(clusters), " clusters X ", length(stim_lab),
                 " stimulation types = ", total_combinations,"."))}

  # Initialize empty output data frames.
  # df_stim_out <- data.frame(matrix(ncol = length(cols), nrow = 0))
  # df_stim_out <- data.frame()
  df_main_out <- data.frame()
  # df_summary_out <- as_tibble(data.frame(matrix(ncol = 3 + length(state_markers), nrow = 0)))
  stacked_bar_plot_data <- data.frame()
  if(umap){
    umap_plot_data <- data.frame()
  }
  df_all_f_out <- data.frame()
  df_all_k_out <- data.frame()
  # k_clust_data <- data.frame()

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
      # df_out_temp <- data.frame(matrix(ncol = length(cols), nrow = 0))
      df_out_temp <- data.frame()

      # Filter cells for a stim type - cell population combinations + unstimulated cells.
      dat_stim_unstim <- dat[dat$stim_type == stim | dat$stim_type %in% unstim_lab, ]
      dat_stim_unstim <- dat_stim_unstim[dat_stim_unstim[,cluster_col] == cluster,]

      # Identify the index of the stimulated and unstimulated cells in the filtered data frame.
      stim_idx <- dat_stim_unstim$stim_type == stim
      unstim_idx <- dat_stim_unstim$stim_type %in% unstim_lab

      # From dat_stim_unstim data frame select expression values for all the state markers
      # and carry out k-means (k = 2) clustering.
      if(verbose == TRUE){message(paste("Carrying out k-means clustering on cells from", cluster, "-", stim, "+ unstim."))}
      dat_state <- dat_stim_unstim[, state_markers]
      k_results <- kmeans(dat_state, 2, iter.max = 100, nstart = 2)
      if (k_results$ifault==4){
        if(verbose == TRUE){message("Using MacQueen algorithm.")}
        k_results <- kmeans(dat_state, 2, iter.max = 100, nstart = 2, algorithm="MacQueen")
      }

      # Record k-means summary for all the combinations.
      df_all_k_out <- rbind(df_all_k_out, cbind(data.frame("stim_type" = stim, "cell_population" = cluster),
                                                broom::tidy(k_results)))

      # Get cluster IDs for stim and unstim cells.
      clust_stim <- k_results$cluster[stim_idx]
      clust_unstim <- k_results$cluster[unstim_idx]

      # Fisher's Exact Test.
      if(verbose == TRUE){message("Carrying out Fisher's exact test.")}

      # Create a contingency table and carry out Fisher's exact test.
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
      con_tab[is.na(con_tab)] = 0
      f_test <- fisher.test(con_tab) # Fisher's exact test.
      f_p_val <- f_test$p.value

      # Record Fisher's exact test statistics for all the tests.
      f_stats <- broom::tidy(f_test)
      df_all_f_out <- rbind(df_all_f_out, cbind(data.frame("stim_type" = stim, "cell_population" = cluster,
                                                           "stim_clust1" = con_tab[1,1], "stim_clust2" = con_tab[2,1],
                                                           "unstim_clust1" = con_tab[1,2], "unstim_clust2" = con_tab[2,2]),
                                                f_stats))

      # Select cells and save data only for the combinations that pass the Fisher's exact test.
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
        # stim_cells_resp <- dat_stim_unstim[k_results$cluster == stim_clust,]
        # stim_cells_resp <- stim_cells_resp[stim_cells_resp$stim_type == stim, ]

        # Combine responding stimulated cells in the output data frame.
        # df_stim_out <- rbind(df_stim_out, stim_cells_resp)

        # Generate data frame with statmarkers and k-means cluster identity per cell.
        f_comb_no <- f_comb_no + 1
        df_main_temp <- dat_stim_unstim %>%
          dplyr::mutate("k_cluster_id" = k_results$cluster) %>%
          dplyr::rename("cell_population" = all_of(cluster_col)) %>%
          dplyr::mutate("responding_cluster" = as.integer(stim_clust)) %>%
          mutate("response_status" = case_when(stim_type == stim & k_cluster_id == responding_cluster ~ "Resp. Stim.",
                                               stim_type == stim & k_cluster_id != responding_cluster ~ "Non-resp. Stim.",
                                               stim_type %in% unstim_lab & k_cluster_id == responding_cluster ~ "Resp. Unstim.",
                                               stim_type %in% unstim_lab & k_cluster_id != responding_cluster ~ "Non-resp. Unstim.")) %>%
          dplyr::mutate("comb_no" = as.integer(f_comb_no))
        df_main_out <- rbind(df_main_temp, df_main_out)

        # Generate data for the summary table.
        # Select median marker expression for each state marker for the selected stim cells.
        # tm <- stim_cells_resp[stim_cells_resp$stim_type == stim, state_markers]
        # med_exp <- apply(tm, 2, median)
        # med_exp <-  as.data.frame(t(med_exp))

        # Calculate fold change from non-responsive to responsive state for stim cells.
        # if(stim_clust == 1){
        #   fc <- round((stim_1/stim_2), digits = 2)
        # }else if(stim_clust == 2){
        #   fc <- round((stim_2/stim_1), digits = 2)
        # }

        # Generate a summary table.
        # summary_table <- cbind(cluster = cluster, stim_type = stim, f_p_value = f_p_val, fold_change = fc, med_exp)
        # df_summary_out <- rbind(df_summary_out, summary_table)

        # Stacked + percent bar plots with stim and unstim on the x-axis and cell counts
        # in k-means clusters on the y-axis.
        stacked_data <- data.frame(stim_status = c("unstim", "unstim", "stim", "stim"),
                                   k_cluster = c("cluster1", "cluster2", "cluster1", "cluster2"),
                                   count = c(unstim_1_, unstim_2_, stim_1_, stim_2_))

        # stacked_temp <- cbind("cell_population" = cluster, "stim_type" = stim, "f_p_val" = f_p_val,
        #                       "stim_clust" = stim_clust, "fold_change" = fc, stacked_data)
        stacked_temp <- cbind("cell_population" = cluster, "stim_type" = stim, "f_p_val" = f_p_val,
                              "stim_clust" = stim_clust, stacked_data)
        stacked_bar_plot_data <- rbind(stacked_bar_plot_data, stacked_temp)

        # Generate UMAP on combined stim and unstim cells and colour them
        # accordingly.
        if(umap == TRUE){
          # Select stimulated cells that are clustered in un-stimulated cells cluster.
          # These are the cells from the simulated samples that likely didn't respond
          # to the stimulant.
          no_of_cells <- min(umap_cells, nrow(df_main_temp))
          set.seed(seed_val)
          input_umap <-  df_main_temp %>% dplyr::slice_sample(n = no_of_cells, replace = FALSE)

          tot_of_cells <- nrow(df_main_temp)

          # Run UMAP.
          if(verbose == TRUE){message(paste0("Running UMAP for ", cluster, " - ", stim,"."))}
          u_res <- uwot::umap(dplyr::select(input_umap, all_of(state_markers)))

          umap_temp <- cbind("cell_population" = cluster, "stim_type" = stim,
                             "tot_of_cells" = tot_of_cells, "no_of_cells" = no_of_cells,
                             "UMAP1" = u_res[,1], "UMAP2" = u_res[,2],
                             "response_status" = as.character(input_umap$response_status))
          umap_plot_data <- rbind(umap_plot_data, umap_temp)
          umap_plot_data$UMAP1 <- as.numeric(umap_plot_data$UMAP1)
          umap_plot_data$UMAP2 <- as.numeric(umap_plot_data$UMAP2)
        }

      } else {
        if(verbose == TRUE){message("Fisher' exact test non-sgnificant. Skipping further steps.")}
      }
    }
  }

  # Adjust p-values for the Fisher's exact test using all tests.
  df_all_f_out <- df_all_f_out %>% mutate("f_p_adj" = p.adjust(p.value, method = "BH"))

  # Generate return list.
  # return_list <- list("selected_expr_data" = df_main_out, "summary" = as_tibble(df_summary_out),
  #                     "stacked_bar_plot_data" = as_tibble(stacked_bar_plot_data),
  #                     "state_markers" = state_markers, "cluster_col" = cluster_col,
  #                     "stim_label" = stim_lab, "unstim_label" = unstim_lab,
  #                     "seed_val" = seed_val, "fisher_test_fail" = as_tibble(df_f_fail_out),
  #                     "all_fisher_p_val" = df_all_f_out, "k_clust_data" = k_clust_data)
  return_list <- list("response_mapping_main" = df_main_out, "stacked_bar_plot_data" = as_tibble(stacked_bar_plot_data),
                      "state_markers" = state_markers, "cluster_col" = cluster_col,
                      "stim_label" = stim_lab, "unstim_label" = unstim_lab,
                      "seed_val" = seed_val, "all_fisher_p_val" = df_all_f_out,
                      "all_k_means_data" = df_all_k_out)
  if(umap == TRUE){
    return_list[["umap_plot_data"]] <- as_tibble(umap_plot_data)
    return_list[["umap"]] <- umap
    return_list[["umap_cells"]] <- umap_cells
  }

 if(verbose == TRUE){message(paste("Response mapping process finished on", date()))}
 return(return_list)
}


