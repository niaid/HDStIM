#' @title HDStIM: High Dimensional Stimulation Immune Mapping
#' @description Function to select cells from the stimulated samples that have likely responded to the stimulant.
#'
#' @param dat              A tibble with the single cell data. Cells on rows
#'                         and variables/markers on columns.
#' @param state_markers    A character vector with the labels of state
#'                         markers from the stimulation panel.
#' @param cellpop_col      Column in the tibble with the cell population IDs.
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
#' @import dplyr
#' @importFrom stats fisher.test kmeans median p.adjust
#' @importFrom uwot umap
#' @importFrom parallel detectCores
#' @export
#'
#' @examples
#' \dontrun{
#' mapped_data <- HDStIM(chi11$expr_data, chi11$state_markers,
#'                   chi11$cluster_col, chi11$stim_label,
#'                   chi11$unstim_label, seed_val = 123, umap = FALSE, umap_cells = NULL,
#'                   verbose = FALSE)
#'}
HDStIM <- function(dat, state_markers, cellpop_col, stim_lab, unstim_lab,
                   seed_val = NULL, umap = FALSE, umap_cells = NULL,
                   verbose = FALSE){
  # For debugging.
  #library(HDStIM)
  # library(tidyverse)
  # dat <- chi11$expr_data
  # state_markers <- chi11$state_markers
  # cellpop_col <- chi11$cluster_col
  # stim_lab <- chi11$stim_label
  # unstim_lab <- chi11$unstim_label
  # seed_val <- 123
  # umap <- TRUE
  # umap_cells <- 50
  # verbose <- TRUE

  # Check argument accuracy.
  if(umap == TRUE & is.null(umap_cells)){
    stop("If umap is set to TRUE then please pass a value for umap_cells.")
  }

  if(verbose == TRUE){message(paste("Mapping started on", date()))}

  # Bind global variables
  p.value <- NULL

  # Standardize column names and state marker labels.
  if(verbose == TRUE){message("Relacing - with _ in state marker labels and expression data column names.")}
  state_markers <- gsub("-", "_", state_markers)
  cols <- colnames(dat)
  cols <- gsub("-", "_", cols)
  colnames(dat) <- cols

  # Convert input data to a tibble.
  dat <- as_tibble(dat)

  # Fetch cluster labels from the expression data frame.
  cellpops <- as.character(unique(dat[[cellpop_col]]))

  total_combinations <- length(cellpops) * length(stim_lab)
  if(verbose == TRUE){message(paste0("Total combinations to go through: ", length(cellpops), " cell populations X ", length(stim_lab),
                 " stimulation types = ", total_combinations,"."))}

  # Initialize empty output data frames.
  df_main_out <- data.frame()
  stacked_bar_plot_data <- data.frame()
  if(umap){
    umap_plot_data <- data.frame()
  }
  df_all_f_out <- data.frame()
  df_all_k_out <- data.frame()

  # Set counter for the number of combinations processed.
  counter <- 0
  f_comb_no <- 0

  # Main nested for loop.
  for(stim in stim_lab){
    for(cellpop in cellpops){
      # Set seed value for the k-means clustering.
      set.seed(seed_val)
      counter <- counter + 1
      if(verbose == TRUE){message(paste0("## Combination ", counter,"/",total_combinations,": ", cellpop, " - ", stim))}
      # Initialize an empty temporary output data frame.
      df_out_temp <- data.frame()

      # Filter cells for a stim type - cell population combinations + unstimulated cells.
      dat_stim_unstim <- dat[dat$stim_type == stim | dat$stim_type %in% unstim_lab, ]
      dat_stim_unstim <- dat_stim_unstim[dat_stim_unstim[,cellpop_col] == cellpop,]

      # Identify the index of the stimulated and unstimulated cells in the filtered data frame.
      stim_idx <- dat_stim_unstim$stim_type == stim
      unstim_idx <- dat_stim_unstim$stim_type %in% unstim_lab

      # From dat_stim_unstim data frame select expression values for all the state markers
      # and carry out k-means (k = 2) clustering.
      if(verbose == TRUE){message(paste("Carrying out k-means clustering on cells from", cellpop, "-", stim, "+ unstim."))}
      dat_state <- dat_stim_unstim[, state_markers]
      k_results <- kmeans(dat_state, 2, iter.max = 100, nstart = 2)
      if (k_results$ifault==4){
        if(verbose == TRUE){message("Using MacQueen algorithm.")}
        k_results <- kmeans(dat_state, 2, iter.max = 100, nstart = 2, algorithm="MacQueen")
      }

      # Record k-means summary for all the combinations.
      df_all_k_out <- rbind(df_all_k_out, cbind(data.frame("stim_type" = stim, "cell_population" = cellpop),
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

      stim_1 <- (as.numeric(length(clust_stim[clust_stim == 1])) / total_stim_count) * 100
      #stim_1_ <- stim_1/total_stim_count
      stim_2 <- (as.numeric(length(clust_stim[clust_stim == 2])) / total_stim_count) * 100
      #stim_2_ <- stim_2/total_stim_count

      unstim_1 <- (as.numeric(length(clust_unstim[clust_unstim == 1])) / total_unstim_count) * 100
      #unstim_1_ <- unstim_1/total_unstim_count
      unstim_2 <- (as.numeric(length(clust_unstim[clust_unstim == 2])) / total_unstim_count) * 100
      #unstim_2_ <- unstim_2/total_unstim_count

      #con_tab <-  matrix(c(round(stim_1_ * 100), round(stim_2_ * 100), round(unstim_1_ * 100), round(unstim_2_ * 100)), nrow = 2, ncol = 2, dimnames = list(c("Cluster1", "Cluster2"), c("Stim", "Unstim")))
      con_tab <-  matrix(c(round(stim_1), round(stim_2), round(unstim_1), round(unstim_2)),
                         nrow = 2, ncol = 2, dimnames = list(c("Cluster1", "Cluster2"), c("Stim", "Unstim")))
      con_tab[is.na(con_tab)] = 0
      f_test <- fisher.test(con_tab) # Fisher's exact test.
      f_p_val <- f_test$p.value

      # Record Fisher's exact test statistics for all the tests.
      f_stats <- broom::tidy(f_test)
      df_all_f_out <- rbind(df_all_f_out, cbind(data.frame("stim_type" = stim, "cell_population" = cellpop,
                                                           "stim_clust1" = con_tab[1,1], "stim_clust2" = con_tab[2,1],
                                                           "unstim_clust1" = con_tab[1,2], "unstim_clust2" = con_tab[2,2]),
                                                f_stats))

      # Select cells and save data only for the combinations that pass the Fisher's exact test.
      if(f_p_val < 0.05){
        if(verbose == TRUE){message("Fisher's exact test significant.")}

        # Identify responding stimulated cells cluster.
        clust_1_ratio <- stim_1 / unstim_1
        clust_2_ratio <- stim_2 / unstim_2
        if(clust_1_ratio > clust_2_ratio){
            unstim_clust = 2
            stim_clust = 1
          } else {
            unstim_clust = 1
            stim_clust = 2
        }

        # Generate data frame with statmarkers and k-means cluster identity per cell.
        f_comb_no <- f_comb_no + 1
        df_main_temp <- dat_stim_unstim %>%
          dplyr::mutate("k_cluster_id" = k_results$cluster) %>%
          dplyr::rename("cell_population" = all_of(cellpop_col)) %>%
          dplyr::mutate("responding_cluster" = as.integer(stim_clust)) %>%
          mutate("response_status" = case_when(stim_type == stim & k_cluster_id == responding_cluster ~ "Resp. Stim.",
                                               stim_type == stim & k_cluster_id != responding_cluster ~ "Non-resp. Stim.",
                                               stim_type %in% unstim_lab & k_cluster_id == responding_cluster ~ "Resp. Unstim.",
                                               stim_type %in% unstim_lab & k_cluster_id != responding_cluster ~ "Non-resp. Unstim.")) %>%
          dplyr::mutate("comb_no" = as.integer(f_comb_no))
        df_main_out <- rbind(df_main_temp, df_main_out)

        # Stacked + percent bar plots with stim and unstim on the x-axis and cell counts
        # in k-means clusters on the y-axis.
        stacked_data <- data.frame(stim_status = c("unstim", "unstim", "stim", "stim"),
                                   k_cluster = c("cluster1", "cluster2", "cluster1", "cluster2"),
                                   cell_count_perc = c(unstim_1, unstim_2, stim_1, stim_2))

        stacked_temp <- cbind("cell_population" = cellpop, "stim_type" = stim, "f_p_val" = f_p_val,
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
          if(verbose == TRUE){message(paste0("Running UMAP for ", cellpop, " - ", stim,"."))}

          n_threads <- parallel::detectCores() - 1
          u_res <- uwot::umap(dplyr::select(input_umap, all_of(state_markers)), n_threads = n_threads)

          umap_temp <- cbind("cell_population" = cellpop, "stim_type" = stim, "condition" = as.character(input_umap$condition),
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
  return_list <- list("response_mapping_main" = df_main_out, "stacked_bar_plot_data" = as_tibble(stacked_bar_plot_data),
                      "state_markers" = state_markers, "cellpop_col" = cellpop_col,
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


