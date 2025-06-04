#' Initialize ToDcall Class with Omics Data
#'
#' This function initializes an object of the ToDcall class using omics data
#' provided directly in the form of data frames. Both the transcript and proteome
#' data frames must have the first column containing gene IDs to allow for
#' matching later. Make sure that they are no duplicates in that column for the transcript count data. Additionally, other column names in both data frames must
#' start with the defined timepoints, followed by a separator (e.g., "_") and
#' a replicate identifier. For example, if timepoints are defined as
#' c("0h", "1h"), then the column names should be like "0h_rep1", "0h_rep2", etc.
#'
#' @param dataset_name A character string for the dataset name.
#' @param organism A character string for the organism name.
#' @param transcript_data A data frame containing the raw transcript counts data and ensembl ids as the first row called "ensembl_base".
#' @param proteome_data A data frame containing the proteome counts data.
#' @param transcript_meta_data A data frame containing the transcript metadata.
#' @param timepoints A vector of timepoints to compare.
#' @param parameters A vector of parameters for ToDcall.
#' @param proteome_normalized A boolean indicating if proteome data is normalized.
#' @return An object of class ToDcall.
#' @export
initToDcall <- function(dataset_name, organism, transcript_data, proteome_data, transcript_meta_data, timepoints, proteome_normalized = T) {

  # Function to check if timepoints match column names
  check_timepoints <- function(data, timepoints) {
    # Extract column names excluding the first column (assumed to be gene IDs)
    col_names <- colnames(data)[-1]

    # Create a pattern to match timepoints at the start of column names with or without replicates
    pattern_with_replicate <- paste0("^(", paste(timepoints, collapse = "|"), ")_")
    pattern_exact <- paste0("^(", paste(timepoints, collapse = "|"), ")$")

    # Check if all columns match either pattern
    all(sapply(col_names, function(name) grepl(pattern_with_replicate, name) || grepl(pattern_exact, name)))
  }

  # Check timepoints for transcript and proteome data
  if (!check_timepoints(transcript_data, timepoints)) {
    stop("Column names in transcript_data do not match the provided timepoints or expected format.")
  }

  if (!check_timepoints(proteome_data, timepoints)) {
    stop("Column names in proteome_data do not match the provided timepoints or expected format.")
  }


  if (proteome_normalized) {
    proteome_norm <- proteome_data
    proteome_counts_raw <- data.frame() # Placeholder, as data is normalized
  } else {
    proteome_counts_raw <- proteome_data
    proteome_norm <- data.frame() # Placeholder, as data is raw
  }

  # Create the ToDcall object
  tod_call_object <- new("ToDcall",
                         dataset_name = dataset_name,
                         organism = organism,
                         transcript_counts_raw = transcript_data,
                         proteome_counts_raw = proteome_counts_raw,
                         transcript_meta_data = transcript_meta_data,
                         timepoints_to_compare = timepoints,
                         proteome_norm = proteome_norm
  )

  return(tod_call_object)
}
