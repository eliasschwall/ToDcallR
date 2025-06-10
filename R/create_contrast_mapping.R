#' Create Contrast Mapping for Time Points
#'
#' This function generates a list of contrasts comparing consecutive time points
#' from a given vector of time points. Each contrast is represented as a vector
#' with the format \code{c("timepoint", tp2, tp1)}, where \code{tp2} and \code{tp1}
#' are consecutive time points.
#'
#' @param timepoints A character vector specifying the time points to compare (e.g., \code{c("0h", "1h", "6h", "12h")}).
#'                   These should be ordered in the sequence they occur.
#' @return A named list where each element is a contrast comparing two consecutive time points.
#'         The names of the list elements are in the format \code{"tp2_vs_tp1"}.
#' @examples
#' timepoints <- c("0h", "1h", "6h", "12h", "24h", "36h", "48h", "72h")
#' contrasts <- create_contrast_mapping(timepoints)
#' print(contrasts)
create_contrast_mapping <- function(timepoints) {
  # Initialize an empty list to store contrasts
  contrasts <- list()

  # Iterate over the timepoints to create contrasts
  for (i in 1:(length(timepoints) - 1)) {
    # Define the current and next timepoints
    tp1 <- timepoints[i]
    tp2 <- timepoints[i + 1]

    # Create a contrast name
    contrast_name <- paste(tp2, "vs", tp1, sep = "_")

    # Add the contrast to the list
    contrasts[[contrast_name]] <- c("timepoint", tp2, tp1)
  }

  return(contrasts)
}
