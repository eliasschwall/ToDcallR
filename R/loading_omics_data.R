loading_omics_data <- function(
    path_to_transcript_data,
    path_to_proteome_data
    ) {

  transcripts <- read.csv(path_to_transcript_data)
  proteome <- read.csv(path_to_proteome_data)

  omics_layers <- list(transcripts, proteome)

  return(omics_layers)
}
