#' @title Marker Ranking ~ Boruta
#' @description Function to run Boruta on the stim - cell population combinations
#'              that passed the Fisher's exact test to rank the markers
#'              according to their contribution to the stim cell selection.
#' @param selected_data      Selected data from the stim_cell_selector function.
#' @param path               Path to the folder for the figures.
#' @param max_runs           Maximum runber of runs for the random forest algorithm.
#'                           Default is 100.
#' @param seed_val           Seed value for Boruta. Default is 123.
#' @param verbose            0, 1, or 2. Default is 0.
#' @return                   Figures saved in a given path.
#'                           And returns a tibble with attributes statistics calculated
#'                           by Boruta.
#' @import Boruta ggplot2
#' @importFrom tibble as_tibble
#' @importFrom grDevices png dev.off rainbow
#' @importFrom graphics plot
#' @export
#'
#' @examples
#' \dontrun{
#' attribute_stats <- marker_ranking_boruta(selected_data, path = NULL,
#'                                         max_runs = 1000, seed_val = 123,
#'                                         verbose = FALSE)
#' }
marker_ranking_boruta <- function(selected_data, path = NULL, max_runs = 100, seed_val = 123, verbose = 0){

  # For debugging.
  # library(Boruta)
  # library(tidyverse)
  # max_runs <- 100
  # verbose <- 2
  # seed_val <- 123
  # path <- file.path("/Users/farmerr2/sandbox/devel/niaid/figures/boruta")

  # Bind global variables.
  comb_no <- maxImp <- medianImp <- minImp <- normHits <- state_marker <- NULL

  state_markers <- selected_data$state_markers
  dat_boruta <- selected_data$k_clust_data
  form_boruta <- as.formula(paste0("k_cluster_id ~ ",paste0(state_markers, collapse = " + ")))

  grouped <- group_by(dat_boruta, comb_no)
  split_data <- group_split(grouped)

  df_stats_out <- data.frame(matrix(nrow = 0, ncol = 9))

  for(i in 1:length(split_data)){
    # dat <- split_data[[1]] # For debugging.
    dat <- split_data[[i]]

    set.seed(seed_val)
    res_boruta <- Boruta(form_boruta, data = dat, doTrace = verbose, maxRuns = max_runs)
    att_stats <- attStats(res_boruta)
    att_stats <- rownames_to_column(att_stats, var = "state_marker")

    stim <- as.character(unique(dat$stim_type)[unique(dat$stim_type) %in% selected_data$stim_label])
    clust <- as.character(unique(dat$cluster))
    att_stats <- cbind("stim_type" = stim, "cluster" = clust, att_stats)
    df_stats_out <- rbind(df_stats_out, att_stats)

    # Plot attribute plot.
    file_name = file.path(path, paste0("imp_", stim, "_", clust, ".png"))
    plot_title = paste0("Stim Type: ", stim, "; Cluster: ", clust)

    png(filename = file_name, res = 300, width = 7, height = 5, units = "in")
    plot(res_boruta, las = 2, cex.axis = 0.7, cex.lab = 0.8, main = plot_title)
    dev.off()

    # Plot attributes using ggplot.
    plot_data <- att_stats[order(att_stats$medianImp),]
    plot_data$state_marker <- factor(plot_data$state_marker, levels = plot_data$state_marker)
    att_plot <- ggplot(plot_data, aes(x = state_marker, y = medianImp, col = normHits)) +
      geom_point() +
      geom_errorbar(aes(ymin=minImp, ymax=maxImp), width = 0.2) +
      labs(x = "State Marker", y = "Importance", col = "Selection Hits",
           title = plot_title) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_colour_gradientn(colours=rainbow(7), values = c(0,0.25, 0.5,0.75,1))

    att_file <- paste0("rank_",stim, "_", clust, ".png")
    ggsave(att_file, plot = att_plot, path = path,
           device = "png", dpi = 300, width = 7, height = 6, units = "in")
  }

  return(as_tibble(df_stats_out))
}
