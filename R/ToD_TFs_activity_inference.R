ToD_TFs_activity_inference <- function(ToDcall) {

  if (ToDcall@organism %in% c("Mouse", "mouse", "Mus Musculus", "mus musculus", "MM", "mm", "Mm")) {
    TFs <- read.delim("tests/testthat/testdata/Mus_musculus_TF.txt")
    organism <- "mouse"
  } else if (ToDcall@organism %in% c("Human", "human", "Homo sapiens", "homo sapiens", "HS", "hs", "Hs")) {
    TFs <- read.delim("tests/testthat/testdata/Homo_sapiens_TF.txt")
    organism <- "human"
  } else {
    stop("Invalid organism specified. Please specify 'Mouse' or 'Human'.")
  }

  ToD_TFs <- lapply(ToDcall@ToD_candidates_filtered, function(time_interval){
    ToD_TFs <- TFs %>% dplyr::filter(Symbol %in% time_interval$mgi_symbol)
    return(ToD_TFs$Symbol)
  }) %>% unlist()

  TFs_and_targets_network <- decoupleR::get_collectri(organism = organism)


  # counts have ensembl ids but need symbols for decoupleR
  ensembl <- biomaRt::useMart("ensembl")
  mouse <- biomaRt::useDataset("mmusculus_gene_ensembl", mart = ensembl)

  # Query to get MGI symbols and Entrez IDs
  transcirptomic_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                                         filters = "ensembl_gene_id",
                                         values = ToDcall@transcript_counts_norm %>% rownames(),
                                         mart = mouse)

  counts <- ToDcall@transcript_counts_norm %>%
    tibble::rownames_to_column("ensembl_gene_id") %>%
    dplyr::left_join(transcirptomic_genes, by = "ensembl_gene_id") %>%
    dplyr::select(-ensembl_gene_id) %>%
    dplyr::mutate(row_mean = rowMeans(select(., -mgi_symbol))) %>%
    dplyr::group_by(mgi_symbol) %>%
    dplyr::filter(row_mean == max(row_mean)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-row_mean) %>%
    dplyr::filter(!is.na(mgi_symbol)) %>%
    tibble::column_to_rownames("mgi_symbol")

  sample_acts <- decoupleR::run_ulm(
    mat= counts,
    net=TFs_and_targets_network,
    .source='source',
    .target='target',
    .mor='mor',
    ) %>%
    dplyr::filter(source %in% ToD_TFs)


  # Transform to wide matrix
  sample_acts_mat <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'condition', names_from = 'source', values_from = 'score') %>%
    tibble::column_to_rownames('condition') %>%
    as.matrix()

  # Get top tfs with more variable means across clusters
  tfs <- sample_acts %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(std = sd(score)) %>%
    dplyr::arrange(-abs(std)) %>%
    dplyr::slice_head(n = n_tfs) %>%  # Use slice_head for head(n_tfs)
    dplyr::pull(source)

  sample_acts_mat <- sample_acts_mat %>% scale() # Scale per sample

  # Choose color palette
  palette_length <- 100
  my_color <- grDevices::colorRampPalette(c("Darkblue", "white", "red"))(palette_length)

  my_breaks <- c(seq(-3, 0, length.out = ceiling(palette_length / 2) + 1),
                 seq(0.05, 3, length.out = floor(palette_length / 2)))

  # Remove underscores from row names
  rownames(sample_acts_mat) <- gsub("_", " ", rownames(sample_acts_mat))

  # Plot
  pheatmap::pheatmap(
    sample_acts_mat,
    border_color = NA,
    color = my_color,
    angle_col = 0,
    main = "ToD TFs activity across samples",
    breaks = my_breaks
  )


}
