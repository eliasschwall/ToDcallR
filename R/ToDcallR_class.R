#' @importClassesFrom DESeq2 DESeqDataSet
setClass("ToDcall", slots = c(
  dataset_name = "character",
  organism = "character",
  transcript_counts_raw = "data.frame",
  proteome_counts_raw = "data.frame",
  transcript_meta_data = "data.frame",
  timepoints_to_compare = "character",
  transcript_counts_norm = "data.frame",
  proteome_norm = "data.frame",
  dds = "DESeqDataSet",
  transcripts_LFC_between_timepoints = "data.frame",
  proteome_LFC_between_timepoints = "data.frame",
  transcripts_LFC_against_0 = "data.frame",
  proteome_LFC_against_0 = "data.frame",
  background_genes_IDs = "data.frame",
  ToDcall_parameters = "numeric",
  ToD_candidates = "list"
))
