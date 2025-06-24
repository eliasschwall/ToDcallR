LFC_Ratio_plot <- function(ToDcall) {

  average_LFC_transcripts <- ToDcall@transcripts_LFC_between_timepoints %>%
    dplyr::group_by(gene_id, comparison) %>% # Group data by gene_id and comparison
    dplyr::summarise(average_LFC = mean(LFC)) %>%
    dplyr::inner_join(ToDcall@background_genes_IDs, by = c("gene_id" = "ensembl_gene_id")) %>%
    dplyr::select(mgi_symbol, comparison, average_LFC) # Select the desired columns

  average_LFC_protein <- ToDcall@proteome_LFC_between_timepoints %>%
    dplyr::group_by(gene_id, comparison) %>% # Group data by gene_id and comparison
    dplyr::summarise(average_LFC = mean(LFC)) # Calculate the average LFC for each group

  # Combine data with left join
  combined_data <- average_LFC_transcripts %>%
    dplyr::left_join(average_LFC_protein, by = c("mgi_symbol" = "gene_id", "comparison")) %>%
    dplyr::rename(
      average_LFC_transcript = average_LFC.x,
      average_LFC_protein = average_LFC.y
    ) %>%
    dplyr::mutate(
      ratio = average_LFC_protein - average_LFC_transcript,
      comparison = factor(comparison, levels = names(ToDcall@ToD_candidates))
    )

  ggplot2::ggplot(combined_data, ggplot2::aes(x = comparison, y = ratio, fill = comparison)) +
    ggplot2::geom_violin() + # Create the violin plot
    ggplot2::ylim(-5,5) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") + # Add cutoff line for increased expression
    ggplot2::geom_hline(yintercept = -1, linetype = "dashed", color = "blue") + # Add cutoff line for decreased expression
    ggplot2::labs(
      title = "Distribution of Protein-RNA LFC-Ratios",
      x = "",
      y = "LFC-Ratio"
    ) +
    ggplot2::theme_minimal() + # Apply a minimal theme
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), # Rotate x-axis labels for better readability
      legend.position = "none"
    ) +
    ggplot2::scale_fill_brewer(palette = "Set3") +
    ggplot2::scale_x_discrete(labels = function(x) {
      x <- gsub("_", " ", x)  # Replace underscores with spaces
      x <- sub("(\\d+h) vs (\\d+h)", "\\2 to \\1", x)  # Swap and replace
      return(x)
    })

}
