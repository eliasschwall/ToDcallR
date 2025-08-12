ToD_TFs_activity_inference <- function(ToDcall) {

  if (ToDcall@organism %in% c("Mouse", "mouse", "Mus Musculus", "mus musculus", "MM", "mm", "Mm")) {
    TFs <- utils::read.delim("tests/testthat/testdata/Mus_musculus_TF.txt")
    organism <- "mouse"
  } else if (ToDcall@organism %in% c("Human", "human", "Homo sapiens", "homo sapiens", "HS", "hs", "Hs")) {
    TFs <- utils::read.delim("tests/testthat/testdata/Homo_sapiens_TF.txt")
    organism <- "human"
  } else {
    stop("Invalid organism specified. Please specify 'Mouse' or 'Human'.")
  }

  ToD_TFs <- lapply(ToDcall@ToD_candidates_filtered, function(time_interval) {
    ToD_TFs <- TFs %>% dplyr::filter(Symbol %in% time_interval$mgi_symbol)
    return(ToD_TFs$Symbol)
  })

  annotation_data <- ToD_TFs %>%
    utils::stack() %>%
    dplyr::rename(
      TFs = values,
      ToD_interval = ind
    ) %>%
    dplyr::mutate(ToD_interval = factor(ToD_interval, levels = unique(ToD_interval)))


  ToD_TFs <- unlist(ToD_TFs)

  TFs_and_targets_network <- decoupleR::get_collectri(organism = organism)

  transcriptomic_genes <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "mgi_symbol"),
    filters = "ensembl_gene_id",
    values = rownames(ToDcall@transcript_counts_norm),
    mart = ToDcall@mart
  )

  counts <- ToDcall@transcript_counts_norm %>%
    tibble::rownames_to_column("ensembl_gene_id") %>%
    dplyr::left_join(transcriptomic_genes, by = "ensembl_gene_id") %>%
    dplyr::select(-ensembl_gene_id) %>%
    dplyr::mutate(row_mean = rowMeans(dplyr::select(., -mgi_symbol))) %>%
    dplyr::group_by(mgi_symbol) %>%
    dplyr::filter(row_mean == max(row_mean)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-row_mean) %>%
    dplyr::filter(!is.na(mgi_symbol)) %>%
    dplyr::distinct(mgi_symbol, .keep_all = T) %>%
    tibble::column_to_rownames("mgi_symbol")

  sample_acts <- decoupleR::run_ulm(
    mat = counts,
    net = TFs_and_targets_network,
    .source = 'source',
    .target = 'target',
    .mor = 'mor'
  ) %>%
    dplyr::filter(source %in% ToD_TFs)

  sample_acts_mat <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'condition', names_from = 'source', values_from = 'score') %>%
    as.data.frame() %>%
    dplyr::mutate(order = as.numeric(sub("h.*", "", condition))) %>%
    dplyr::arrange(order) %>%
    dplyr::select(-order) %>%
    tibble::column_to_rownames('condition') %>%
    as.matrix()


  sample_acts_mat <- scale(sample_acts_mat)

  palette_length <- 100
  my_color <- grDevices::colorRampPalette(c("Darkblue", "white", "red"))(palette_length)

  # Adjust breaks to match the color length
  my_breaks <- seq(-3, 3, length.out = palette_length)

  rownames(sample_acts_mat) <- gsub("_", " ", rownames(sample_acts_mat))

  # filtering annotation data to only include TFs from CollecTRI
  annotation_data <- annotation_data %>%
    dplyr::filter(TFs %in% colnames(sample_acts_mat))

  # Extract unique genes and time intervals
  genes <- unique(annotation_data$TFs)
  time_intervals <- unique(annotation_data$ToD_interval)

  # Create a binary matrix indicating the presence of each gene in each time interval
  presence_matrix <- matrix(0, nrow = length(genes), ncol = length(time_intervals))
  rownames(presence_matrix) <- genes
  colnames(presence_matrix) <- time_intervals

  # Fill the matrix dynamically based on annotation_data
  for (i in seq_len(nrow(annotation_data))) {
    gene <- annotation_data$TFs[i]
    interval <- annotation_data$ToD_interval[i]
    if (gene %in% genes && interval %in% time_intervals) {
      presence_matrix[gene, as.character(interval)] <- 1
    }
  }

  # Get the row names from both matrices
  sample_acts_row_names <- rownames(t(sample_acts_mat))
  presence_matrix_row_names <- rownames(presence_matrix)

  # Reorder presence_matrix according to the order in sample_acts_mat
  presence_matrix <- presence_matrix[match(sample_acts_row_names, presence_matrix_row_names), ]

  colnames(presence_matrix) <- sapply(colnames(presence_matrix), function(x) {
    x <- gsub("_", " ", x)  # Replace underscores with spaces
    x <- sub("(\\d+h) vs (\\d+h)", "\\2 to \\1", x)  # Swap and replace
    return(x)
  })

  annotation_heatmap <- ComplexHeatmap::Heatmap(
    presence_matrix,
    name = "ToD Interval",
    col = c("0" = "white", "1" = "lightgreen"),
    show_row_names = FALSE,
    rect_gp = grid::gpar(col = "black", lwd = 2),
    show_column_names = T,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_rot = 45,
    heatmap_legend_param = list(
      title = "ToD Regulation",
      border = "black",
      at = c(0, 1),  # Specify the values for the legend
      labels = c("No ToD", "ToD")  # Labels for each value
    )
  )

  # Plot using ComplexHeatmap
  main_heatmap <- ComplexHeatmap::Heatmap(
    t(sample_acts_mat),
    name = "Activity",
    rect_gp = grid::gpar(col = "black", lwd = 2),
    col = circlize::colorRamp2(my_breaks, my_color),
    cluster_rows = TRUE,
    row_names_side = "left",
    cluster_columns = FALSE,
    column_title = "ToD TFs activity across samples",
    column_title_gp = grid::gpar(fontsize = 27),
    column_names_rot = 45
  )

  ComplexHeatmap::draw(main_heatmap + annotation_heatmap)

}


