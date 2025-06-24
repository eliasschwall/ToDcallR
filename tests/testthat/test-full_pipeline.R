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

ToDcall <- calculate_LFCs(ToDcall, proteome_log = T, transcriptome_log = F)

ToDcall_strict <- callToDs(ToDcall,stable_transcript_range = 1, ToD_threshold = 1, ToD_filtering_next_time_point_range = 0.3)
ToDcall_loose <- callToDs(ToDcall,stable_transcript_range = 2, ToD_threshold = 1, ToD_filtering_next_time_point_range = 0.3)

enrichments_strict <- ORA(ToDcall_strict, 0.1)
enrichments_loose <- ORA(ToDcall_loose, 0.1)

A <- LFC_Ratio_plot(ToDcall)
B <- ToDCs_barplot(ToDcall, plot_filtered_list = F, ylim = 120)
C <- ToDCs_barplot(ToDcall, plot_filtered_list = T, ylim = 120)

library(patchwork)
patchwork <- A / ((B + C) + plot_layout(guides = 'collect'))
yang_tod_number_plot <- patchwork + plot_annotation(tag_levels = "A")

ggplot2::ggsave(filename = "yang_tod_plot.svg", plot = yang_tod_number_plot, device = "svg")


#rounding padj
enrichments_loose$enrichments$`6h_vs_1h`$MH@result$p.adjust <- round(enrichments_loose$enrichments$`6h_vs_1h`$MH@result$p.adjust,3)

p1 <- enrichments_loose$enrichments$`1h_vs_0h`$MH %>% enrichplot::dotplot(title = "0h to 1h Hallmarks")
p2 <- enrichments_loose$enrichments$`6h_vs_1h`$MH %>% enrichplot::dotplot(title = "1h to 6h Hallmarks")
p3 <- enrichments_loose$enrichments$`6h_vs_1h`$`M5_GO:BP` %>% enrichplot::dotplot(title = "1h to 6h GO Biological Processes")
p4 <- enrichments_loose$enrichments$`12h_vs_6h`$`M5_GO:MF` %>% enrichplot::dotplot(title = "6h to 12h GO Molecular Function")

ORA_plot <- (p1+p2)/(p3+p4) + plot_annotation(tag_levels = "A")
ggplot2::ggsave(filename = "/Users/elias/Documents/ComputationalBiologyMaster/Master_Thesis/figures/ORA_plot.svg", plot = ORA_plot, device = "svg", width = 15, height = 10)



TFs_mouse <- read.delim("tests/testthat/testdata/Mus_musculus_TF.txt")
