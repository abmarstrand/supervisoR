#' Load and preprocess gene sets
#'
#' @importFrom stringr str_to_title str_remove str_replace_all
#' @importFrom igraph graph_from_data_frame V vertex edges
#' @importFrom dplyr select distinct mutate rename
#' @importFrom msigdbr msigdbr
#' @importFrom readr read_tsv
#' @importFrom magrittr %>%
#' @param database Which database should we load data for? Currently supports 'reactome' and 'GO' corresponding to the msigdb Reactome or GOBP databases.
#' @param species Which species should we load data for? Support species in the msigdb database, currently 'Homo sapiens' and 'Mus musculus'.
#' @return A list containing gene sets, mapping data frame, and the pathway graph.
#' @export
load_default_data <- function(database = "reactome",
                              species = "Homo sapiens") {
  if (!database %in% c("reactome", "GO")) {
    stop(paste0("Database: ", database, " not supported by default."))
  }
  # Set parameters based on the database
  params <- list(
    reactome = list(
      path = system.file("extdata", "ReactomePathwaysRelation.txt", package = "supervisoR"),
      prefix = "^REACTOME ",
      category = "C2",
      subcategory = "CP:REACTOME"
    ),
    GO = list(
      path = system.file("extdata", "GOTermRelation.txt", package = "supervisoR"),
      prefix = "^GOBP ",
      category = "C5",
      subcategory = "GO:BP"
    )
  )
  db_params <- params[[database]]

  # Load msigdbr data
  msigdb_init <- msigdbr(
    species = species,
    category = db_params$category,
    subcategory = db_params$subcategory) %>%
    select(gs_name, gs_exact_source, gene_symbol) %>%
    distinct() %>%
    rename(name = gs_name, exact_source = gs_exact_source)

  processed_names <- msigdb_init$name %>%
    str_replace_all("_", " ") %>%
    str_remove(db_params$prefix) %>%
    str_to_title()

  msigdb_init <- msigdb_init %>%
    mutate(processed_name = processed_names)

  # Create geneset
  geneset <- split(msigdb_init$gene_symbol, msigdb_init$processed_name)

  # Create mapping
  mapping <- msigdb_init %>%
    select(processed_name, exact_source) %>%
    distinct()

  # Read the parent-child relationships
  relations <- read_tsv(db_params$path, col_names = TRUE, show_col_types = FALSE)

  # Ensure 'relation' column is present
  # If not present (Reactome), create it and set to NA
  if (!"relation" %in% names(relations)) {
    relations$relation <- NA_character_
  }

  # Identify intermediate nodes: appear in relations but not in mapping$exact_source
  all_nodes <- unique(c(relations$parent, relations$child))
  known_nodes <- mapping$exact_source
  intermediates <- setdiff(all_nodes, known_nodes)

  # If there are intermediates, flatten the hierarchy in a vectorized manner
  if (length(intermediates) > 0) {
    # Identify edges where child is intermediate
    child_is_intermediate <- relations$child %in% intermediates
    # Identify edges where parent is intermediate
    parent_is_intermediate <- relations$parent %in% intermediates
    
    # Split parents by intermediate child
    # parents_by_intermediate[[node]] = all parents of that intermediate node
    parents_by_intermediate <- split(relations$parent[child_is_intermediate],
                                     relations$child[child_is_intermediate])
    parents_by_intermediate <- lapply(parents_by_intermediate, unique)
    
    # Split children by intermediate parent
    # children_by_intermediate[[node]] = all children of that intermediate node
    children_by_intermediate <- split(relations$child[parent_is_intermediate],
                                      relations$parent[parent_is_intermediate])
    children_by_intermediate <- lapply(children_by_intermediate, unique)
    
    # Remove all edges involving intermediate nodes
    keep_edges <- !(relations$parent %in% intermediates | relations$child %in% intermediates)
    relations <- relations[keep_edges, , drop = FALSE]
    
    # Create a list to store new edges
    node_names <- intersect(names(parents_by_intermediate), names(children_by_intermediate))
    total_expansions <- 0L
    np_vec <- integer(length(node_names))
    nc_vec <- integer(length(node_names))
    
    for (i in seq_along(node_names)) {
      node_parents <- parents_by_intermediate[[ node_names[i] ]]
      node_children <- children_by_intermediate[[ node_names[i] ]]
      if (!is.null(node_children) && length(node_parents) > 0 && length(node_children) > 0) {
        np <- length(node_parents)
        nc <- length(node_children)
        expansions <- np * nc
        np_vec[i] <- np
        nc_vec[i] <- nc
        total_expansions <- total_expansions + expansions
      } else {
        np_vec[i] <- 0L
        nc_vec[i] <- 0L
      }
    }
    
    if (total_expansions > 0) {
      # Pre-allocate the vectors
      all_parents <- character(total_expansions)
      all_children <- character(total_expansions)
      all_relations <- rep(NA_character_, total_expansions)
      
      # Fill them
      pos <- 1L
      for (i in seq_along(node_names)) {
        np <- np_vec[i]
        nc <- nc_vec[i]
        if (np > 0 && nc > 0) {
          node_parents <- parents_by_intermediate[[ node_names[i] ]]
          node_children <- children_by_intermediate[[ node_names[i] ]]
          expansions <- np * nc
          
          # Create parent/child sequences
          # parent repeated for each child's group
          # child repeated for each parent's group
          parents_block <- rep(node_parents, times = nc)
          children_block <- rep(node_children, each = np)
          
          # Fill the allocated space
          all_parents[pos:(pos + expansions - 1L)] <- parents_block
          all_children[pos:(pos + expansions - 1L)] <- children_block
          # all_relations is already NA_character_
          pos <- pos + expansions
        }
      }
      
      # Create a single data frame of new edges
      new_edges <- data.frame(
        parent = all_parents,
        child = all_children,
        relation = all_relations,
        stringsAsFactors = FALSE
      )
      
      # Combine and ensure uniqueness
      relations <- unique(rbind(relations, new_edges))
    }
  }
  # Build the final graph
  g <- graph_from_data_frame(relations, directed = TRUE)
  # Map exact_source IDs to processed names
  exact_source_to_name <- setNames(mapping$processed_name, mapping$exact_source)
  # Assign labels to graph vertices
  V(g)$label <- exact_source_to_name[V(g)$name]
  # Assign genes to vertices
  V(g)$genes <- geneset[V(g)$label]
  V(g)$genes[sapply(V(g)$genes, is.null)] <- list(character(0))
  return(list(
    gene_sets = geneset,
    mapping = mapping,
    pathway_graph = g,
    relations = relations
  ))
}

#' Load and preprocess gene sets and pathway relations
#'
#' @importFrom igraph graph_from_data_frame V vertex edges
#' @importFrom dplyr left_join rename
#' @importFrom stats setNames
#' @importFrom stringr str_to_title
#' @param gene_sets Named list of gene sets, with pathways as names and genes as elements.
#' @param pathways_relation Data frame with two columns: 'parent' and 'child' pathways.
#' @param translation_layer Optional data frame with two columns: 'gene_sets_name' and 'relation_name'.
#' @return A list containing gene sets, mapping data frame and the pathway graph.
#' @export
load_and_preprocess_gene_sets <- function(gene_sets,
                                          pathways_relation,
                                          translation_layer = NULL) {
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
    mapping <- left_join(mapping,
                         translation_layer,
                         by = c("processed_name" = "gene_sets_name"))
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

  all_exact_sources <- unique(c(pathways_relation$parent,
                                pathways_relation$child))
  complete_mapping <- data.frame(
    exact_source = all_exact_sources,
    stringsAsFactors = FALSE
  )
  complete_mapping <- left_join(complete_mapping,
                                mapping,
                                by = "exact_source")

  complete_mapping$processed_name <- ifelse(
    is.na(complete_mapping$processed_name),
    ifelse(grepl("REACTOME", complete_mapping$exact_source[1]),
           str_to_title(gsub("_", " ", gsub("^REACTOME_", "", complete_mapping$exact_source))),
           str_to_title(gsub("_", " ", gsub("^GOBP_", "", complete_mapping$exact_source)))
    ),
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
  V(g)$genes <- lapply(V(g)$genes,
                       function(x) if (is.null(x)) character(0) else x)

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

  from_genes_list <- V(g)$genes[match(from_nodes, V(g)$name)]
  to_genes_list <- V(g)$genes[match(to_nodes, V(g)$name)]

  # Compute overlaps and percent overlaps
  overlaps <- mapply(function(from, to) length(intersect(from, to)),
                     from_genes_list, to_genes_list)
  from_sizes <- lengths(from_genes_list)
  percent_overlaps <- ifelse(from_sizes > 0, (overlaps / from_sizes) * 100, 0)

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
  # Check that parent_geneset_name is length 1 and valid
  if (length(parent_geneset_name) == 0) {
    stop("No parent_geneset_name provided or it is an empty vector.")
  }

  # Get the exact_source ID for the parent pathway
  parent_indices <- match(parent_geneset_name, mapping$processed_name)

  # If no match was found
  if (length(parent_indices) == 0 || is.na(parent_indices)) {
    stop("Parent pathway not found in the mapping.")
  }

  parent_exact_source <- mapping$exact_source[parent_indices]

  # Check if parent_exact_source is missing
  if (length(parent_exact_source) == 0 || is.na(parent_exact_source)) {
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
  descendant_geneset_names <- mapping$processed_name[match(descendant_exact_sources, mapping$exact_source)]
  descendant_geneset_names <- descendant_geneset_names[!is.na(descendant_geneset_names)]

  return(descendant_geneset_names)
}
