get_msigdb_gene_sets <- function(organism) {

  if (organism %in% c("Mouse", "mouse", "Mus Musculus", "mus musculus", "MM", "mm", "Mm")) {
    spec <- "Mus musculus"
    db <- "MM"
  } else if (organism %in% c("Human", "human", "Homo sapiens", "homo sapiens", "HS", "hs", "Hs")) {
    spec <- "Homo sapiens"
    db <- "HS"
  } else {
    stop("Invalid organism specified. Please specify 'Mouse' or 'Human'.")
  }

  collection_mapping <- msigdbr::msigdbr_collections(db) %>%
    dplyr::group_by(gs_collection) %>%
    dplyr::summarize(subcollections = list(unique(gs_subcollection))) %>%
    tibble::deframe()

  genesets <- lapply(names(collection_mapping), function(collection) {
    subcollections <- collection_mapping[[collection]]

    # Iterate over each subcollection for the current collection
    subcollection_results <- lapply(subcollections, function(subcollection) {
      # Determine the subcollection argument
      subcollection_arg <- if (subcollection == "") NULL else subcollection

      # Retrieve the data frame for the given collection and subcollection
      msigdbr::msigdbr(
        species = spec,
        db_species = db,
        collection = collection,
        subcollection = subcollection_arg
      )
    })

    # Name each subcollection result
    names(subcollection_results) <- ifelse(subcollections == "", collection, paste(collection, subcollections, sep = "_"))

    return(subcollection_results)
  })

  # Flatten the list of results if needed
  genesets <- do.call(c, genesets)

}
