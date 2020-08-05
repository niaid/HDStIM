#' @title State (signalling) Marker Selection
#' @description Function to select state (signalling) markers that best discriminate
#'              stimulated cells from unstimulated cells per stimulation condition.
#'              This method is based on a combination of linear mixed effects modelling
#'              and ANOVA. The selected list of markers can be used in the multi-marker
#'              stim cell selection function.
#' @param dat              A tibble with the single cell data. Cells on rows
#'                         and variables/markers on columns.
#' @param state_markers    A character vector with the labels of state
#'                         markers from the stimulation panel.
#' @param stim_lab         A character vector of stim label(s).
#' @param unstim_lab       A character of unstim label(s).
#' @param REML             Logical. To use restricted maximum likely hood estimation
#'                         in model fitting.
#' @param verbose          Logical. To make function more verbose. Default is FALSE.
#'
#' @return A tibble with p-value from ANOVA per stim type and state marker.
#' @importFrom tibble as_tibble
#' @import dplyr
#' @importFrom stats anova
#' @importFrom lme4 lmer
#' @export
#'
#' @examples
#' selected_state_markers <- state_marker_selection(chi11_1k$expr_data, chi11_1k$state_markers,
#'                  chi11_1k$cluster_col, chi11_1k$stim_label, chi11_1k$unstim_label,
#'                  REML = FALSE, verbose = FALSE)
state_marker_selection <- function(dat, state_markers, cluster_col, stim_lab, unstim_lab, REML = FALSE, verbose = FALSE){
  #for debugging
  # dat <- chi11_1k$expr_data
  # state_markers <- chi11_1k$state_markers
  # cluster_col <- chi11_1k$cluster_col
  # stim_lab <- chi11_1k$stim_label
  # unstim_lab <- chi11_1k$unstim_label
  # REML <- FALSE
  # verbose <- TRUE

  df_out <- data.frame(matrix(nrow = 0, ncol = 3))
  for(stim in stim_lab){
    dat_stim_unstim <- dat[dat$stim_type == stim | dat$stim_type == unstim_lab,]
    for(state in state_markers){
      if(verbose == TRUE){message(paste("Computing lmer for", stim, "&", state))}
      form <- as.formula(paste(state, " ~ ", "stim_type", " + (1 | ", cluster_col, ")"))
      lmm <- lmer(form, data = dat_stim_unstim, REML = REML)

      form_ran <- as.formula(paste(state, " ~ (1 | ", cluster_col, ")"))
      lmm_ran <- lmer(form_ran, data = dat_stim_unstim, REML = REML)

      if(verbose == TRUE){message(paste("Computing ANOVA for", stim, "&", state))}
      if(REML){anv <- anova(lmm_ran, lmm, refit = FALSE)}
      else{anv <- anova(lmm_ran, lmm)}
      df_out <- rbind(df_out, data.frame("stim_type" = stim, "state_marker" = state, "anova_p_val" = anv$Pr[2]))
    }
  }

  # Do FDR correction.
  grouped <- group_by(df_out, stim_type)
  df_out_adj <- mutate(grouped, anova_p_adj = p.adjust(anova_p_val, method='fdr'))
  return(ungroup(df_out_adj))
}

