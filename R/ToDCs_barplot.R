ToDCs_barplot <- function(ToDcall) {

  ToD_candidates_per_time_interval <- ToDcall@ToD_candidates_filtered %>%
    lapply(nrow) %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(cols = everything(), names_to = "comparison", values_to = "count") %>%
    dplyr::mutate(
      comparison = factor(comparison, levels = names(ToDcall@ToD_candidates)),
      end_time = as.numeric(stringr::str_extract(comparison, "^[0-9]+")),
      start_time = as.numeric(stringr::str_extract(comparison, "(?<=_vs_)[0-9]+")),
      time_interval = end_time - start_time,
      normalized_count = count / time_interval
      ) %>%
    dplyr::select(comparison, count, normalized_count) %>%
    tidyr::pivot_longer(cols = c(count, normalized_count), names_to = "type", values_to = "value")


  ggplot2::ggplot(ToD_candidates_per_time_interval, ggplot2::aes(x = comparison, y = value, fill = type)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(
      title = "Number of ToD candidates per time interval",
      x = "",
      y = "ToD Candidates",
      fill = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = c(0.05, 0.99), # Position in the top left inside the plot
      legend.justification = c("left", "top") # Align the legend
    ) +
    ggplot2::scale_fill_manual(
      values = c("count" = "darkgreen", "normalized_count" = "lightgreen"),
      labels = c("count" = "Count", "normalized_count" = "Normalized for time") # Custom labels
      ) +
    ggplot2::scale_x_discrete(labels = function(x) gsub("_", " ", x))

}
