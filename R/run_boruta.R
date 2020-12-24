#' @title Run Boruta
#' @description Function to run Boruta on the stim - cell population combinations
#'              that passed the Fisher's exact test to rank the markers
#'              according to their contribution to the stim cell selection.
#' @param selected_data
#'
#' @return
#' @export
#'
#' @examples
run_boruta <- function(selected_data, verbose = 0){

  state_markers <- selected_data$state_markers
  dat_boruta <- selected_data$k_clust_data
  form_boruta <- as.formula(paste0("k_cluster_id ~ ",paste0(state_markers, collapse = " + ")))

  grouped <- group_by(dat_boruta, stim_type, cluster)
  split_data <- group_split(grouped)

  df_stats_out <- data.frame(matrix(nrow = 0, ncol = 9))
  for(dat in split_data){
    # dat <- split_data[[1]]
    res_boruta <- Boruta(form_boruta, data = dat, doTrace = verbose, maxRuns = 1000)
    att_stats <- attStats(res_boruta)
    att_stats <- rownames_to_column(att_stats, var = "state_marker")
    att_stats <- cbind("stim_type" = unique(dat$stim_type), "cluster" = unique(dat$cluster), att_stats)
    df_stats_out <- rbind(df_stats_out, att_stats)
  }

}
plot_boruta <- plot(res_boruta)
