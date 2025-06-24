UTR_analysis <- function(ToDcall) {

  # fetching the UTRs for ToD candidates
  ToDC_UTR_metrics <- lapply(names(ToDcall@ToD_candidates_filtered), function(time_interval) {
    fetch_UTRs(ToDcall@ToD_candidates_filtered[[time_interval]], mart.object = ToDcall@mart)
  })

  # set the names of the list to the time intervals
  names(ToDC_UTR_metrics) <- names(ToDcall@ToD_candidates_filtered)

  # fetching the UTRs for the background genes
  UTRs_metrics_background_genes <- fetch_UTRs(ToDcall@background_genes_IDs, ToDcall@mart)

  # convert them to Biostrings
  background_5_UTRs_biostring <- Biostrings::DNAStringSet(UTRs_metrics_background_genes$utr_5$five_prime_utr)
  names(background_5_UTRs_biostring) <- UTRs_metrics_background_genes$utr_5_metrics$mgi_symbol

  background_3_UTRs_biostring <- Biostrings::DNAStringSet(UTRs_metrics_background_genes$utr_3$three_prime_utr)
  names(background_3_UTRs_biostring) <- UTRs_metrics_background_genes$utr_3_metrics$mgi_symbol

  # converting ToDC_UTR to Biostrings as well
  ToDC_5_UTRs_biostring <- lapply(ToDC_UTR_metrics, function(interval) {

    utr_5_biostring <- interval$utr_5_metrics$five_prime_utr %>%
      Biostrings::DNAStringSet() %>%
      setNames(., interval$utr_5_metrics$mgi_symbol)

    return(utr_5_biostring)
  })

  ToDC_3_UTRs_biostring <- lapply(ToDC_UTR_metrics, function(interval) {

    utr_3_biostring <- interval$utr_3_metrics$three_prime_utr %>%
      Biostrings::DNAStringSet() %>%
      setNames(., interval$utr_3_metrics$mgi_symbol)

    return(utr_3_biostring)
  })

  ToDcall@ToD_UTRs_analysis <- list(
    ToDC_UTR_metrics = ToDC_UTR_metrics,
    ToDC_5_UTRs_biostring = ToDC_5_UTRs_biostring,
    ToDC_3_UTRs_biostring = ToDC_3_UTRs_biostring,
    UTRs_metrics_background_genes = UTRs_metrics_background_genes,
    background_5_UTRs_biostring = background_5_UTRs_biostring,
    background_3_UTRs_biostring = background_3_UTRs_biostring
  )

  return(ToDcall)
}
