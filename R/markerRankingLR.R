#' @title Marker Ranking ~ Logistic Regression
#' @description            Logistic regression for each marker to test it's contribution in
#'                         k-means clustering.
#' @param selected_data    Selected data from the stim_cell_selector function.
#' @param lr_max_it        Maximum iteration (default = 50) for logistic regression.
#' @param verbose          Boolean (TRUE/FALSE). Default is FALSE.
#' @importFrom stats glm coefficients
#' @export
#'
#' @examples
#' \dontrun{
#' lr_stats <- marker_ranking_lr(selected_data, lr_max_it = 50, verbose = FALSE)
#' }
marker_ranking_lr <- function(selected_data, lr_max_it = 50, verbose = FALSE){

  # For Debugging
  # lr_max_it <- 50
  # verbose <- TRUE

  # Bind global variables.
  comb_no <- state_marker <- NULL

  state_markers <- selected_data$state_markers
  dat_for_lr <- selected_data$k_clust_data

  df_lr_out <- data.frame(matrix(nrow = 0, ncol = 8))

  grouped <- group_by(dat_for_lr, comb_no)
  split_data <- group_split(grouped)
  for(i in 1:length(split_data)){
    dat_for_lr <- split_data[[i]]
    stim <- as.character(unique(dat_for_lr$stim_type)[unique(dat_for_lr$stim_type) %in% selected_data$stim_label])
    cluster <- as.character(unique(dat_for_lr$cluster))
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

 return(as_tibble(df_lr_out))
}

