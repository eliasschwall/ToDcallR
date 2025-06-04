test_that("initToDcall works", {

  library(dplyr)
  sample_names <- colnames(read_csv("tests/testthat/testdata/transcriptome.csv"))[-1]
  colData <- dplyr::tibble(sample = sample_names) %>%
    dplyr::mutate(
      # Extract timepoint (e.g., "X0h" becomes "0h")
      timepoint = str_extract(sample, "\\d+h"),
      # Extract replicate info (e.g., "rep1")
      replicate = str_extract(sample, "rep\\d+")
    ) %>%
    # Convert to factors with a specified level order for timepoint.
    dplyr::mutate(
      timepoint = factor(timepoint, levels = c("0h", "1h", "6h", "12h", "24h", "36h", "48h", "72h")),
      replicate = factor(replicate)
    ) %>%
    tibble::column_to_rownames("sample")


  ToDcall <- initToDcall(
    dataset_name = "Yang Data",
    organism = "Mouse",
    transcript_data = read_csv("tests/testthat/testdata/transcriptome.csv"),
    proteome_data = read_csv("tests/testthat/testdata/proteome.csv"),
    proteome_normalized = T,
    timepoints = c("0h","1h","6h","12h","24h","36h","48h","72h"),
    transcript_meta_data = colData
  )

  expect_s4_class(ToDcall, "ToDcall")

})
