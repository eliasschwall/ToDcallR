callToDs <- function(ToDcall, stable_transcript_range, ToD_threshold, ToD_filtering_next_time_point_range) {

  ToDcall@ToDcall_parameters <- c(stable_transcript_range, ToD_threshold, ToD_filtering_next_time_point_range) %>% setNames(c("stable_transcript_range", "ToD_threshold","ToD_filtering_next_time_point_range"))

  # get gene ids for background genes
  if (ToDcall@organism %in% c("Mouse", "mouse", "Mus Musculus", "mus musculus", "MM", "mm", "Mm")) {
    mart <- biomaRt::useDataset("mmusculus_gene_ensembl", mart = biomaRt::useMart("ensembl"))
  } else if (ToDcall@organism %in% c("Human", "human", "Homo sapiens", "homo sapiens", "HS", "hs", "Hs")) {
    mart <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = biomaRt::useMart("ensembl"))
  } else {
    stop("Invalid organism specified. Please specify 'Mouse' or 'Human'.")
  }

  ToDcall@mart <- mart

  # Query to get MGI symbols and Entrez IDs
  transcirptomic_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "entrezgene_id"),
                            filters = "ensembl_gene_id",
                            values = ToDcall@transcript_counts_norm %>% rownames(),
                            mart = mart)

  SummarizedExperiment::rowData(ToDcall@dds)$mgi_symbol <- transcirptomic_genes$mgi_symbol[match(rownames(ToDcall@dds), transcirptomic_genes$ensembl_gene_id)]


  ToDcall@background_genes_IDs <- transcirptomic_genes %>%
    dplyr::filter(mgi_symbol %in% ToDcall@proteome_norm$gene_id)

  # define contrast to extract stable genes
  contrasts <- create_contrast_mapping(ToDcall@timepoints_to_compare)

  # filter for stable genes
  stable_genes <- list()

  for (comparison in names(contrasts)) {
    # Run DESeq2 results with a modified null hypothesis:
    # Testing if the true log₂FC is within ±1 (null hypothesis) versus outside that range.
    stable_genes[[comparison]] <- DESeq2::results(
      ToDcall@dds,
      contrast = contrasts[[comparison]],
      lfcThreshold = stable_transcript_range,
      altHypothesis = "lessAbs"
      ) %>%
      as.data.frame() %>%
      dplyr::filter(padj<0.05) %>%
      tibble::rownames_to_column("ensembl_gene_id") %>%
      dplyr::left_join(transcirptomic_genes, by = "ensembl_gene_id")

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
      dplyr::filter(diff >= ToD_threshold) %>%
      dplyr::filter(!is.na(mgi_symbol))

    return(ToD_candidates)
  })

  names(ToD_candidates) <- names(stable_genes)


  ToD_candidates_filtered <- lapply(names(ToD_candidates), function(comp) {
    # Get the index of the current comparison
    current_index <- match(comp, names(ToD_candidates))

    # If it's the last comparison, we cannot check the next time point
    if (current_index == length(names(ToD_candidates))) {
      return(NULL)
    }

    # Get the next comparison
    next_comp <- names(ToD_candidates)[current_index + 1]

    # Determine the current and next comparison names in ToDcall@proteome_LFC_against_0
    current_comp_mapped <- sub("h_vs_.*", "h_vs_0h", comp)
    next_comp_mapped <- sub("h_vs_.*", "h_vs_0h", next_comp)

    # Filter the current candidates
    ToD_candidates[[comp]] %>%
      dplyr::filter(mgi_symbol %in% ToDcall@proteome_LFC_against_0$gene_id) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        current_LFC = ToDcall@proteome_LFC_against_0 %>%
          dplyr::filter(gene_id == mgi_symbol, comparison == current_comp_mapped) %>%
          dplyr::pull(LFC) %>% dplyr::first(),
        next_LFC = ToDcall@proteome_LFC_against_0 %>%
          dplyr::filter(gene_id == mgi_symbol, comparison == next_comp_mapped) %>%
          dplyr::pull(LFC) %>% dplyr::first()
      ) %>%
      dplyr::filter(next_LFC >= (1-ToD_filtering_next_time_point_range) * current_LFC)  # Allow up to 30% decrease
  })
  names(ToD_candidates_filtered) <- names(stable_genes)


  ToDcall@ToD_candidates <- ToD_candidates
  ToDcall@ToD_candidates_filtered <- ToD_candidates_filtered %>% Filter(Negate(is.null), .)





  return(ToDcall)
}
