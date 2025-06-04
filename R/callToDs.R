callToDs <- function(ToDcall, stable_transcript_range, ToD_threshold) {

  ToDcall@ToDcall_parameters <- c(stable_transcript_range, ToD_threshold) %>% setNames(c("stable_transcript_range", "ToD_threshold"))

  # get gene ids for background genes
  ensembl <- biomaRt::useMart("ensembl")
  mouse <- biomaRt::useDataset("mmusculus_gene_ensembl", mart = ensembl)

  # Query to get MGI symbols and Entrez IDs
  ToDcall@background_genes_IDs <- biomaRt::getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "entrezgene_id"),
                            filters = "ensembl_gene_id",
                            values = ToDcall@transcript_counts_norm %>% rownames(),
                            mart = mouse)

  # define contrast to extract stable genes
  contrasts <- create_contrast_mapping(ToDcall@timepoints_to_compare)

  # filter for stable genes
  stable_genes <- list()

  for (comparison in names(contrasts)) {
    # Run DESeq2 results with a modified null hypothesis:
    # Testing if the true log₂FC is within ±1 (null hypothesis) versus outside that range.
    res_stable <- DESeq2::results(
      ToDcall@dds,
      contrast = contrasts[[comparison]],
      lfcThreshold = 1,
      altHypothesis = "greaterAbs"
      ) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("ensembl_gene_id") %>%
      dplyr::left_join(ToDcall@background_genes_IDs, by = "ensembl_gene_id")

    # Filter for genes whose 95% confidence interval lies entirely within [stable_transcript_range]
    stable_genes[[comparison]] <- res_stable %>%
      dplyr::filter((log2FoldChange - 1.96 * lfcSE) >= -stable_transcript_range,
                    (log2FoldChange + 1.96 * lfcSE) <= stable_transcript_range)

    cat("Comparison:", comparison, "\n")

  }

  # now we match the stable genes with the proteome LFCs to calculate the LFC ratios and call the ToD candidates
  ToD_candidates <- lapply(names(stable_genes), function(comp){
    stable_transcripts <- stable_genes[[comp]] %>% dplyr::select(mgi_symbol, log2FoldChange)
    protein_LFC <- ToDcall@proteome_LFC_between_timepoints %>%
      dplyr::filter(comparison == comp) %>%
      dplyr::select(gene_id, LFC)

    # Merge the data frames based on mgi_symbol and gene_id and filter for ToD candidates
    ToD_candidates <- merge(stable_transcripts, protein_LFC, by.x = "mgi_symbol", by.y = "gene_id", all = TRUE) %>%
      dplyr::mutate(diff = LFC - log2FoldChange) %>%
      dplyr::filter(diff >= ToD_threshold)

    return(ToD_candidates)
  })

  names(ToD_candidates) <- names(stable_genes)
  ToDcall@ToD_candidates <- ToD_candidates

  return(ToDcall)
}
