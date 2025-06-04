#' Calculate Log Fold Changes (LFC) Compared to the First Time Point
#'
#' This function calculates the log fold changes for each gene compared to the first
#' specified time point, using expression data with or without replicates.
#'
#' @param data A data frame containing gene expression data. The first column should be gene IDs,
#'             and subsequent columns should be expression values. Column names can be formatted
#'             as "time_replicate" (e.g., "0h_rep1") or just "time" (e.g., "0h").
#' @param timepoints A character vector specifying the time points to compare (e.g., c("0h", "1h", "6h", "12h")).
#' @param is_log_transformed A boolean indicating if the input data is already log-transformed.
#' @return A data frame containing the LFCs for each time point compared to the first time point.
#'         The data frame includes columns: gene_id, replicate, comparison, and LFC.
#' @importFrom dplyr inner_join mutate select filter bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @importFrom stats log2
LFC_calc_against_first_timepoint <- function(data, timepoints, is_log_transformed = FALSE) {

  # Check if column names include replicates
  if (any(grepl("_", colnames(data)[-1]))) {
    # If replicates are present, separate time and replicate
    data_long <- tidyr::pivot_longer(data, -gene_id, names_to = c("time", "replicate"), names_sep = "_", values_to = "expression")
  } else {
    # If no replicates, treat the whole column name as time
    data_long <- tidyr::pivot_longer(data, -gene_id, names_to = "time", values_to = "expression")
    data_long$replicate <- "rep1"  # Assign a dummy replicate name
  }

  # Initialize an empty data frame to store results
  lfc_results <- data.frame()

  # Set the first time point
  tp1 <- timepoints[1]

  # Calculate LFC for each time point compared to the first time point
  for (j in 2:length(timepoints)) {
    tp2 <- timepoints[j]

    # Filter data for each time point
    data_tp1 <- dplyr::filter(data_long, time == tp1)
    data_tp2 <- dplyr::filter(data_long, time == tp2)

    # Join data on gene_id and replicate
    lfc_data <- dplyr::inner_join(data_tp1, data_tp2, by = c("gene_id", "replicate"), suffix = c(".tp1", ".tp2"))

    # Calculate LFC
    if (is_log_transformed) {
      lfc_data <- dplyr::mutate(lfc_data, LFC = expression.tp2 - expression.tp1)
    } else {
      lfc_data <- dplyr::mutate(lfc_data, LFC = log2(expression.tp2 + 0.0001) - log2(expression.tp1 + 0.0001))
    }

    # Add a column for the comparison
    lfc_data <- dplyr::mutate(lfc_data, comparison = paste(tp2, "vs", tp1, sep = "_"))

    # Select relevant columns
    lfc_data <- dplyr::select(lfc_data, gene_id, replicate, comparison, LFC)

    # Append to results
    lfc_results <- dplyr::bind_rows(lfc_results, lfc_data)
  }

  return(lfc_results)
}
