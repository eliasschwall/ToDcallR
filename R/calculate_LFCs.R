
#' Title
#'
#' @param ToDcall ToDcallR object
#' @param proteome_log boolean if the proteome data is already log transformed
#' @param transcriptome_log boolean if the transcriptome data is already log transformed
#'
#' @returns
#' @export
#'
#' @examples
#' @importFrom magrittr %>%
calculate_LFCs <- function(ToDcall, proteome_log, transcriptome_log) {

  # First, check if the normalized transcriptomic data is available and if not we normalize the count data
  if (nrow(ToDcall@transcript_counts_norm) == 0) {

    # DESeq normalization
    counts <- ToDcall@transcript_counts_raw %>% tibble::column_to_rownames("ensembl_base")
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      colData = ToDcall@transcript_meta_data,
      design = ~ timepoint
    ) %>% DESeq2::DESeq()
    ToDcall@dds <- dds
    ToDcall@transcript_counts_norm <- DESeq2::counts(dds, normalized = T) %>% as.data.frame()
  }

  # now we do the same for the same for the proteome data
  if(nrow(ToDcall@proteome_norm) == 0){
    # log transform the proteome data
    ToDcall@proteome_norm <- ToDcall@proteome_counts_raw %>% dplyr::mutate(dplyr::across(-1, ~ log2(. + 1)))
  }

  # now we calculate the LFCs for the transcriptome data between timepoints
  ToDcall@transcripts_LFC_between_timepoints <- ToDcall@transcript_counts_norm %>%
    tibble::rownames_to_column("gene_id") %>%
    LFC_calc_between_timepoints(., ToDcall@timepoints_to_compare, is_log_transformed = transcriptome_log)

  # now we calculate the LFCs for the proteome data between timepoints
  ToDcall@proteome_LFC_between_timepoints <- ToDcall@proteome_norm %>%
    dplyr::rename(gene_id = 1) %>%
    LFC_calc_between_timepoints(., ToDcall@timepoints_to_compare, is_log_transformed = proteome_log)

  # now we calculate the LFCs for the transcriptome data against the first timepoint
  ToDcall@transcripts_LFC_against_0 <- ToDcall@transcript_counts_norm %>%
    tibble::rownames_to_column("gene_id") %>%
    LFC_calc_against_first_timepoint(., ToDcall@timepoints_to_compare, is_log_transformed = transcriptome_log)

  # now we calculate the LFCs for the proteome data against the first timepoint
  ToDcall@proteome_LFC_against_0 <- ToDcall@proteome_norm %>%
    dplyr::rename(gene_id = 1) %>%
    LFC_calc_against_first_timepoint(., ToDcall@timepoints_to_compare, is_log_transformed = proteome_log)

  return(ToDcall)
}
