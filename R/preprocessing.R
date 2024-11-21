#' Load and preprocess gene sets
#'
#' @importFrom jsonlite fromJSON
#' @importFrom stringr str_to_title
#' @importFrom igraph graph_from_data_frame V delete_vertices
#' @importFrom magrittr %>%
#' @importFrom dplyr select distinct mutate rename
#' @importFrom msigdbr msigdbr
#' @importFrom readr read_tsv
#' @importFrom purrr map_chr
#' @importFrom rlang .data
#' @return A list containing gene sets, mapping data frame, and the pathway graph.
#' @export
load_default_data <- function() {
  pathways_relation_path <- system.file("extdata",
                                        "ReactomePathwaysRelation.txt",
                                        package = "supervisoR")

  msigdb_init <- msigdbr(species = "Homo sapiens",
                         category = "C2",
                         subcategory = "CP:REACTOME") %>%
    select(.data$gs_name, .data$gs_exact_source, .data$gene_symbol) %>%
    distinct() %>%
    rename(name = .data$gs_name,
           exact_source = .data$gs_exact_source,
           gene_symbol = .data$gene_symbol)

  geneset <- split(msigdb_init$gene_symbol, msigdb_init$name)
  names(geneset) <- str_to_title(gsub("REACTOME ", "", gsub("_", " ", names(geneset))))

  mapping <- msigdb_init %>%
    select(.data$name, .data$exact_source) %>%
    distinct() %>%
    mutate(processed_name = str_to_title(gsub("REACTOME ", "", gsub("_", " ", .data$name)))) %>%
    as.data.frame()

  # Read the parent-child relationships
  relations <- read_tsv(pathways_relation_path, col_names = FALSE)
  colnames(relations) <- c("parent", "child")

  # Build the graph
  g <- graph_from_data_frame(relations, directed = TRUE)

  # Map geneset names to exact_source IDs
  geneset_names <- names(geneset)
  geneset_to_exact_source <- mapping$exact_source[match(geneset_names, mapping$processed_name)]
  names(geneset_to_exact_source) <- geneset_names

  # Create reverse mapping
  exact_source_to_name <- mapping$processed_name
  names(exact_source_to_name) <- mapping$exact_source

  # Assign labels to graph vertices
  V(g)$label <- exact_source_to_name[V(g)$name]
  # Assign default labels to vertices with missing labels
  V(g)$label[is.na(V(g)$label)] <- V(g)$name[is.na(V(g)$label)]
  V(g)$genes <- geneset[exact_source_to_name[V(g)$name]]

  # Assign genes to vertices
  #V(g)$genes <- geneset[exact_source_to_name[V(g)$name]]
  V(g)$genes <- lapply(V(g)$genes, function(x) if (is.null(x)) character(0) else x)

  return(list(
    gene_sets = geneset,
    mapping = mapping,
    pathway_graph = g,
    relations = relations
  ))
}

#' Load and preprocess gene sets and pathway relations
#'
#' @importFrom igraph graph_from_data_frame V
#' @importFrom dplyr left_join rename
#' @importFrom stats setNames
#' @importFrom stringr str_to_title
#' @param gene_sets Named list of gene sets, with pathways as names and genes as elements.
#' @param pathways_relation Data frame with two columns: 'parent' and 'child' pathways.
#' @param translation_layer Optional data frame with two columns: 'gene_sets_name' and 'relation_name'.
#' @return A list containing gene sets, mapping data frame and the pathway graph.
#' @export
load_and_preprocess_gene_sets <- function(gene_sets, pathways_relation, translation_layer = NULL) {
  # Ensure that gene_sets is a named list
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a named list with pathway names as names.")
  }

  # Ensure that pathways_relation is a data frame with 'parent' and 'child' columns
  if (!is.data.frame(pathways_relation) || !all(c("parent", "child") %in% colnames(pathways_relation))) {
    stop("pathways_relation must be a data frame with 'parent' and 'child' columns.")
  }

  # Process the translation layer if provided
  if (!is.null(translation_layer)) {
    if (!is.data.frame(translation_layer) || ncol(translation_layer) != 2) {
      stop("translation_layer must be a data frame with two columns.")
    }
    colnames(translation_layer) <- c("gene_sets_name", "relation_name")
  }

  # Create mapping data frame
  gene_sets_names <- names(gene_sets)
  mapping <- data.frame(
    processed_name = gene_sets_names,
    stringsAsFactors = FALSE
  )

  # If translation_layer is provided, map the names
  if (!is.null(translation_layer)) {
    mapping <- left_join(mapping, translation_layer, by = c("processed_name" = "gene_sets_name"))
    # Use 'relation_name' as the name to match in pathways_relation
    mapping$exact_source <- mapping$relation_name
  } else {
    # Without translation layer, assume names match between gene_sets and pathways_relation
    mapping$exact_source <- mapping$processed_name
  }

  # Handle missing exact_source
  missing_exact_source <- is.na(mapping$exact_source)
  if (any(missing_exact_source)) {
    warning("Some pathways in gene_sets do not have matching entries in the translation_layer. They will be excluded.")
    mapping <- mapping[!missing_exact_source, ]
  }

  all_exact_sources <- unique(c(pathways_relation$parent, pathways_relation$child))
  complete_mapping <- data.frame(
    exact_source = all_exact_sources,
    stringsAsFactors = FALSE
  )
  complete_mapping <- left_join(complete_mapping, mapping, by = "exact_source")

  complete_mapping$processed_name <- ifelse(
    is.na(complete_mapping$processed_name),
    str_to_title(gsub("_", " ", gsub("^REACTOME_", "", complete_mapping$exact_source))),
    complete_mapping$processed_name
  )

  # Create exact_source_to_name mapping
  exact_source_to_name <- setNames(complete_mapping$processed_name, complete_mapping$exact_source)

  # Filter pathways_relation to include only pathways present in the mapping
  valid_exact_sources <- mapping$exact_source
  pathways_relation_filtered <- pathways_relation[pathways_relation$parent %in% valid_exact_sources &
                                                    pathways_relation$child %in% valid_exact_sources, ]

  # Create the graph
  g <- graph_from_data_frame(d = pathways_relation_filtered, directed = TRUE)

  # Assign labels to graph vertices
  V(g)$label <- exact_source_to_name[V(g)$name]
  # Assign default labels to vertices with missing labels
  missing_labels <- is.na(V(g)$label)
  V(g)$label[missing_labels] <- V(g)$name[missing_labels]

  # Assign genes to vertices
  V(g)$genes <- gene_sets[exact_source_to_name[V(g)$name]]
  V(g)$genes <- lapply(V(g)$genes, function(x) if (is.null(x)) character(0) else x)

  return(list(
    gene_sets = gene_sets,
    mapping = mapping,
    pathway_graph = g
  ))
}


#' Add overlap information to edges
#'
#' @importFrom igraph E ends
#' @importFrom dplyr intersect
#' @param g The pathway graph (igraph object).
#' @param gene_sets A named list of gene sets.
#' @return The graph with overlap information added to edges.
#' @export
add_overlap_to_edges <- function(g, gene_sets) {
  # Retrieve edge endpoints
  edge_ends <- ends(g, E(g))
  from_nodes <- edge_ends[, 1]
  to_nodes <- edge_ends[, 2]

  from_genes_list <- gene_sets[V(g)[from_nodes]$label]
  to_genes_list <- gene_sets[V(g)[to_nodes]$label]
  from_genes_list <- lapply(from_genes_list, function(x) if (is.null(x)) character(0) else x)
  to_genes_list <- lapply(to_genes_list, function(x) if (is.null(x)) character(0) else x)


  # Compute overlaps and percent overlaps
  overlaps <- mapply(function(from, to) length(intersect(from, to)), from_genes_list, to_genes_list)
  from_sizes <- sapply(from_genes_list, length)
  percent_overlaps <- (overlaps / from_sizes) * 100
  percent_overlaps[is.na(percent_overlaps)] <- 100
  overlaps[is.na(overlaps)] <- 1

  # Assign edge attributes
  E(g)$overlap <- overlaps
  E(g)$percent_overlap <- percent_overlaps

  return(g)
}

#' Get child pathways of a parent pathway
#'
#' @importFrom igraph subcomponent V
#' @param parent_geneset_name Name of the parent geneset.
#' @param mapping Data frame containing mapping information.
#' @param g The pathway graph (igraph object).
#' @return A character vector of child pathway names.
#' @export
get_child_pathways <- function(parent_geneset_name, mapping, g) {
  # Get the exact_source ID for the parent pathway
  parent_exact_source <- mapping$exact_source[match(parent_geneset_name, mapping$processed_name)]
  if (is.na(parent_exact_source)) {
    stop("Parent pathway not found in the mapping.")
  }

  # Find the parent vertex in the graph
  parent_vertex <- which(V(g)$name == parent_exact_source)
  if (length(parent_vertex) == 0) {
    stop("Parent pathway not found in the graph.")
  }

  # Get all descendants (children) of the parent pathway
  descendants <- subcomponent(g, v = parent_vertex, mode = "out")

  # Exclude the parent itself
  descendant_exact_sources <- V(g)$name[descendants]
  descendant_exact_sources <- descendant_exact_sources[descendant_exact_sources != parent_exact_source]

  # Map exact_source IDs to pathway names
  descendant_geneset_names <- mapping[match(descendant_exact_sources, mapping$exact_source), "processed_name"]
  descendant_geneset_names <- descendant_geneset_names[!is.na(descendant_geneset_names)]

  return(descendant_geneset_names)
}
