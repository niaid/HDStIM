#' @title Sample data set for CyTOF Stimulation Assay
#'
#' @description A \code{list} with the CyTOF stimulation assay data.
#'
#' @format A list with one \code{tibble} containig CyTOF expression data.
#'         And four \code{character vectors} for arguments in the \code{\link{HDStIM}} function.
#' \describe{
#' \item{chi11$expr_data}{A 7,000 X 36 \code{tibble}. Cells are on the rows
#'                          and variables on the columns. The first 6 columns contain
#'                          for each cell \code{cluster_id} (from \code{FlowSOM} clustering),
#'                          \code{sample_id} (unique for each FSC file),
#'                          \code{condition} (comparison groups), \code{patient_id} (unique for each subject),
#'                          \code{stim_type} (labels for types of stimulation assays including the unstim),
#'                          \code{merging1} (meta culster labels from \code{ConsensusClusterPlus}).
#'                          The last 30 columns contain the \code{archsinh} transformed CyTOF expression
#'                          values for the 30 markers (20 type and 10 state) used in the sitmulation panel.}
#' \item{chi11$type_markers}{A character vector with the labels for type markers used in the stimulation panel.}
#' \item{chi11$state_markers}{A character vector with the labels for state markers used in the stimulation panel.}
#' \item{chi11$cluster_col}{A character label of the meta-cluster/cluster ID column in \code{chi11$expr_dat} tibble.}
#' \item{chi11$stim_label}{A character vector with the label(s) for the stimulation types corresponding to the
#'                            labels in the\code{stim_type} column in \code{chi11$expr_data}.}
#' \item{chi11$unstim_label}{A character label for the unstim cells corresponding to the
#'                            labels in the\code{stim_type} column in \code{chi11$expr_data}.}
#' }
"chi11"
