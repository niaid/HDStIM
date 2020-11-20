#' @title State (signalling) Marker Selection Using Linear Mixed Effects Modelling
#' @description Function to select state (signalling) markers that best discriminate
#'              stimulated cells from unstimulated cells per stimulation condition.
#'              This method is based on a combination of linear mixed effects modelling
#'              and ANOVA. The selected list of markers can be used in the multi-marker
#'              stim cell selection function.
#' @param dat              A tibble with the single cell data. Cells on rows
#'                         and variables/markers on columns.
#' @param state_markers    A character vector with the labels of state
#'                         markers from the stimulation panel.
#' @param cluster_col      Name (character) of the column with cluster/cell population names
#'                         in the input tibble.
#' @param stim_lab         A character vector of stim label(s).
#' @param unstim_lab       A character of unstim label(s).
#' @param REML             Logical. To use restricted maximum likely hood estimation
#'                         in model fitting.
#' @param verbose          Logical. To make function more verbose. Default is FALSE.
#'
#' @return A tibble with p-value from ANOVA per stim type and state marker.
#' @importFrom tibble as_tibble rownames_to_column
#' @import dplyr
#' @importFrom stats anova as.formula p.adjust
#' @importFrom lme4 lmer lmerControl .makeCC
#' @export
#'
#' @examples
#' selected_state_markers <- state_marker_selection_lmm(chi11_1k$expr_data, chi11_1k$state_markers,
#'                  chi11_1k$cluster_col, chi11_1k$stim_label, chi11_1k$unstim_label,
#'                  REML = FALSE, verbose = FALSE)
state_marker_selection_lmm <- function(dat, state_markers, cluster_col, stim_lab, unstim_lab, REML = FALSE, verbose = FALSE){
  #for debugging
  # dat <- chi11_1k$expr_data
  # state_markers <- chi11_1k$state_markers
  # cluster_col <- chi11_1k$cluster_col
  # stim_lab <- chi11_1k$stim_label
  # unstim_lab <- chi11_1k$unstim_label
  # REML <- FALSE
  # verbose <- TRUE

  # Bind global variables.
  stim_type <- anova_p_val <- NULL

  df_out <- data.frame(matrix(nrow = 0, ncol = 4))
  df_anova <- data.frame(matrix(nrow = 0, ncol = 15))
  clusters <- unique(dat[[cluster_col]])
  for(stim in stim_lab){
    for(clust in clusters){
      dat_stim_unstim <- dat[dat$stim_type == stim | dat$stim_type == unstim_lab,]
      dat_stim_unstim <- dat_stim_unstim[dat_stim_unstim[[cluster_col]] == clust,]
      for(state in state_markers){
        form <- as.formula(paste(state, " ~ ", "stim_type", " + (1 | sample_id)"))
        lmm <- lmer(form, data = dat_stim_unstim, REML = REML)

        form_ran <- as.formula(paste(state, " ~ (1 | sample_id)"))
        lmm_ran <- lmer(form_ran, data = dat_stim_unstim, REML = REML,
                        control = lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))

        med_stim <- median(dat_stim_unstim[dat_stim_unstim$stim_type == stim, ][[state]]) + 1
        med_unstim <- median(dat_stim_unstim[dat_stim_unstim$stim_type == unstim_lab, ][[state]]) + 1
        fc <- med_stim / med_unstim

        if(verbose == TRUE){message(paste("Computing ANOVA for", stim, "&", state))}
        if(REML){anv <- anova(lmm_ran, lmm, refit = FALSE)}else{anv <- anova(lmm_ran, lmm)}

        df_out <- rbind(df_out, data.frame("stim_type" = stim, "cluster" = clust,
                                           "state_marker" = state, "fc" = fc, "anova_p_val" = anv$Pr[2]))

        # May not be needed in the final version (only for testing)
        no_stim_cells <- nrow(dat_stim_unstim[dat_stim_unstim$stim_type == stim,])
        no_unstim_cells <- nrow(dat_stim_unstim[dat_stim_unstim$stim_type == unstim_lab,])
        df_anv <- rownames_to_column(data.frame(anv), "model")
        df_temp <- cbind("stim_type" = stim, "cluster" = clust, "state_marker" = state,
                         "total_cells" = nrow(dat_stim_unstim), "no_stim_cells" = no_stim_cells,
                         "no_unstim_cells" = no_unstim_cells, df_anv)
        df_anova <- rbind(df_anova, df_temp)
      }
    }
  }

  # Do FDR correction.
  grouped <- group_by(df_out, stim_type)
  df_out_adj <- mutate(grouped, anova_p_adj = p.adjust(anova_p_val, method='fdr'))
  return(list("df_fdr" = ungroup(df_out_adj), "df_anova" = df_anova))
}

# Deprecated Code - DELETE EVENTUALLY
#state_marker_selection <- function(dat, state_markers, cluster_col, stim_lab, unstim_lab, REML = FALSE, verbose = FALSE){
  #for debugging
  # dat <- chi11_1k$expr_data
  # state_markers <- chi11_1k$state_markers
  # cluster_col <- chi11_1k$cluster_col
  # stim_lab <- chi11_1k$stim_label
  # unstim_lab <- chi11_1k$unstim_label
  # REML <- FALSE
  # verbose <- TRUE

#   df_out <- data.frame(matrix(nrow = 0, ncol = 4))
#   clusters <- unique(dat[[cluster_col]])
#   for(stim in stim_lab){
#     for(clust in clusters){
#       dat_stim_unstim <- dat[dat$stim_type == stim | dat$stim_type == unstim_lab,]
#       dat_stim_unstim <- dat_stim_unstim[dat_stim_unstim[[cluster_col]] == clust,]
#       dat_stim_unstim$stim_type <- as.numeric(dat_stim_unstim$stim_type == stim)
#       for(state in state_markers){
#         form <- as.formula(paste0("stim_type", " ~ ", state))
#         fm <- glm(form, data = dat_stim_unstim, family = "binomial")
#         p_val <- coef(summary(fm))[8]
#         df_out <- rbind(df_out, data.frame("stim_type" = stim, "cluster" = clust,  "state_marker" = state, "glm_p_val" = p_val))
#       }
#     }
#   }
#
#   # Do FDR correction.
#   grouped <- group_by(df_out, stim_type)
#   df_out_adj <- mutate(grouped, glm_p_adj = p.adjust(glm_p_val, method='fdr'))
#   return(ungroup(df_out_adj))
# }
