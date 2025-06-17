ORA <- function(ToDcall) {

  # list for all msigdb gene sets that could be relevant
  msigdb_collection <- get_msigdb_gene_sets(ToDcall@organism)

  # run enrichment for each ToD time points for each gene set collection
  enrichments <- lapply(ToDcall@ToD_candidates_filtered, function(time_interval){
    lapply(msigdb_collection, function(genesets){
      clusterProfiler::enricher(
        gene = time_interval$mgi_symbol,
        TERM2GENE = genesets %>% dplyr::select(gs_name, gene_symbol),
        universe = ToDcall@background_genes_IDs$mgi_symbol
      )
    })
  })

  enrichment_plots <- lapply(names(enrichments), function(interval_name) {
    time_interval <- enrichments[[interval_name]]

    # Create a named list for each gene set within the time interval
    gene_set_plots <- lapply(names(time_interval), function(gene_set_name) {
      gene_set_enrichment <- time_interval[[gene_set_name]]

      if (nrow(gene_set_enrichment) > 0) {
        # Create the dot plot with a title including both interval and gene set names
        enrichplot::dotplot(gene_set_enrichment) +
          ggtitle(paste(
            "Dataset:", ToDcall@dataset_name,
            "\nGene Set Collection:", gene_set_name,
            "\nTime Interval:", interval_name,
            "\nStable Transcript Range:", ToDcall@ToDcall_parameters[1],
            "\nToD Threshold:", ToDcall@ToDcall_parameters[2],
            "\nToD filtering next time point range:", ToDcall@ToDcall_parameters[3]
                        ))
      } else {
        NULL  # Explicitly return NULL for clarity, but will be filtered out
      }
    })

    # Assign names to the list of plots for each gene set
    names(gene_set_plots) <- names(time_interval)

    # Filter out NULL plots
    Filter(Negate(is.null), gene_set_plots)
  })

  # Assign names to the list of time intervals
  names(enrichment_plots) <- names(enrichments)

  return(enrichment_plots)
}
