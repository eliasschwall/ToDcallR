ToDCs_barplot <- function(ToDcall, plot_filtered_list = T, ylim = 150) {

  if(plot_filtered_list == T){
    candidates <- ToDcall@ToD_candidates_filtered
    subtitle <- "filtered"
  } else if(plot_filtered_list == F){
    candidates <- ToDcall@ToD_candidates
    subtitle <- "no filtering"
  } else{
    stop("plot_filtered_list must be either True or False")
  }

  ToD_candidates_per_time_interval <- candidates %>%
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
      title = "Number of ToD candidates",
      subtitle = subtitle,
      x = "",
      y = "ToD Candidates",
      fill = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::scale_fill_manual(
      values = c("count" = "darkgreen", "normalized_count" = "lightgreen"),
      labels = c("count" = "Count", "normalized_count" = "Normalized for time") # Custom labels
      ) +
    ggplot2::scale_x_discrete(labels = function(x) {
      x <- gsub("_", " ", x)  # Replace underscores with spaces
      x <- sub("(\\d+h) vs (\\d+h)", "\\2 to \\1", x)  # Swap and replace
      return(x)
    }) +
    ggplot2::ylim(0,ylim)

}
