#' Stim Cell Selector
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
#' @param seed_val         Seed value for k-means clustering.
#'
#' @return
#' @export
#'
#' @examples
#' dat <- as_tibble(read.table(file.path("path/data.tsv"), sep = "\t", headers = T))
#' state_markers <- c("pSTAT1", "pSTAT5", "IkBa", "pCREB")
#' cluster_col <- "cluster_id"
#' stim_lab <- c("A", "T", "L", "G") # A for interferon alpha,
#'                                   # T for TCR, L for LPS and
#'                                   # G for interferon gamma.
#' unstim_lab <- "U"
#' seed_val <- 123
#' selected_cells <- stim_cell_selector(dat, state_markers, cluster_col, stim_lab, unstim_lab, seed_val = seed_val
stim_cell_selector <- function(dat, state_markers, cluster_col, stim_lab, unstim_lab, seed_val = F){

}
