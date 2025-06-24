fetch_UTRs <- function(gene_df, mart.object = mart) {

  # fetching the biomaRt attributes we need separately for both UTRs
  utr_5 <- biomaRt::getBM(attributes = c("5utr", "mgi_symbol", "ensembl_gene_id"),
                          filters = "mgi_symbol",
                          values = gene_df$mgi_symbol,
                          mart = mart.object)

  utr_3 <- biomaRt::getBM(attributes = c("3utr", "mgi_symbol", "ensembl_gene_id"),
                          filters = "mgi_symbol",
                          values = gene_df$mgi_symbol,
                          mart = mart.object)

  # Change column names
  colnames(utr_5)[1] <- "five_prime_utr"
  colnames(utr_3)[1] <- "three_prime_utr"

  # Process 5' UTRs
  utr_5 <- utr_5 %>%
    dplyr::tibble() %>%
    dplyr::filter(!stringr::str_detect(five_prime_utr, "Sequence unavailable")) %>%
    dplyr::mutate(
      utr_5_length = nchar(five_prime_utr),
      A_count_absolute = stringr::str_count(five_prime_utr, "A"),
      A_count_percentage = stringr::str_count(five_prime_utr, "A") / utr_5_length * 100,
      T_count_absolute = stringr::str_count(five_prime_utr, "T"),
      T_count_percentage = stringr::str_count(five_prime_utr, "T") / utr_5_length * 100,
      G_count_absolute = stringr::str_count(five_prime_utr, "G"),
      G_count_percentage = stringr::str_count(five_prime_utr, "G") / utr_5_length * 100,
      C_count_absolute = stringr::str_count(five_prime_utr, "C"),
      C_count_percentage = stringr::str_count(five_prime_utr, "C") / utr_5_length * 100,
      A_to_T_ratio = stringr::str_count(five_prime_utr, "A") / stringr::str_count(five_prime_utr, "T"),
      G_to_C_ratio = stringr::str_count(five_prime_utr, "G") / stringr::str_count(five_prime_utr, "C"),
      GC_content_percentage = ((stringr::str_count(five_prime_utr, "G") + stringr::str_count(five_prime_utr, "C")) /
                                 (stringr::str_count(five_prime_utr, "A") + stringr::str_count(five_prime_utr, "T") +
                                    stringr::str_count(five_prime_utr, "G") + stringr::str_count(five_prime_utr, "C"))*100)
    )

  # Process 3' UTRs
  utr_3 <- utr_3 %>%
    dplyr::tibble() %>%
    dplyr::filter(!stringr::str_detect(three_prime_utr, "Sequence unavailable")) %>%
    dplyr::mutate(
      utr_3_length = nchar(three_prime_utr),
      A_count_absolute = stringr::str_count(three_prime_utr, "A"),
      A_count_percentage = stringr::str_count(three_prime_utr, "A") / utr_3_length * 100,
      T_count_absolute = stringr::str_count(three_prime_utr, "T"),
      T_count_percentage = stringr::str_count(three_prime_utr, "T") / utr_3_length * 100,
      G_count_absolute = stringr::str_count(three_prime_utr, "G"),
      G_count_percentage = stringr::str_count(three_prime_utr, "G") / utr_3_length * 100,
      C_count_absolute = stringr::str_count(three_prime_utr, "C"),
      C_count_percentage = stringr::str_count(three_prime_utr, "C") / utr_3_length * 100,
      A_to_T_ratio = stringr::str_count(three_prime_utr, "A") / stringr::str_count(three_prime_utr, "T"),
      G_to_C_ratio = stringr::str_count(three_prime_utr, "G") / stringr::str_count(three_prime_utr, "C"),
      GC_content_percentage = ((stringr::str_count(three_prime_utr, "G") + stringr::str_count(three_prime_utr, "C")) /
                                 (stringr::str_count(three_prime_utr, "A") + stringr::str_count(three_prime_utr, "T") +
                                    stringr::str_count(three_prime_utr, "G") + stringr::str_count(three_prime_utr, "C"))*100)
    )

  utr_metrics <- list(utr_5_metrics = utr_5, utr_3_metrics = utr_3)

  return(utr_metrics)
}
