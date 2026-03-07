#' Load built-in pathway resources
#'
#' Load the package's built-in pathway gene sets, mapping table, and pathway
#' relations for one supported database/species combination. This is the easiest
#' entry point for pathway analyses and visualization
#'
#' @details
#' The returned object contains everything needed by [plot_subgraph()] and
#' [run_pathway_shiny_app()]. `database = "reactome"` loads Reactome relations
#' bundled with the package, while `database = "GO"` loads Gene Ontology
#' biological process relations. Gene memberships are obtained via `msigdbr` and
#' pathway relations are read from files stored in `inst/extdata/`.
#'
#' @importFrom stringr str_to_title str_remove str_replace_all
#' @importFrom igraph graph_from_data_frame V vertex edges V<-
#' @importFrom dplyr select distinct mutate rename
#' @importFrom msigdbr msigdbr
#' @importFrom readr read_tsv
#' @importFrom magrittr %>%
#' @param database Character scalar identifying the built-in database. Supported
#'   values are `"reactome"` and `"GO"`.
#' @param species Character scalar passed to `msigdbr::msigdbr()`, for example
#'   `"Homo sapiens"` or `"Mus musculus"`.
#' @return A list with four elements: `gene_sets` (named list), `mapping`
#'   (data frame with processed and exact pathway identifiers), `pathway_graph`
#'   (`igraph` object), and `relations` (edge table).
#' @examples
#' \donttest{
#' data <- load_default_data(database = "reactome", species = "Homo sapiens")
#' names(data)
#' head(data$mapping)
#' }
#' @seealso [load_and_preprocess_gene_sets()], [plot_subgraph()]
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
      collection = "C2",
      subcollection = "CP:REACTOME"
    ),
    GO = list(
      path = system.file("extdata", "GOTermRelation.txt", package = "supervisoR"),
      prefix = "^GOBP ",
      collection = "C5",
      subcollection = "GO:BP"
    )
  )
  db_params <- params[[database]]
  
  # msigdbr has version differences - try both.
  msigdb_raw <- tryCatch(
    {
      msigdbr(
        species = species,
        collection = db_params$collection,
        subcollection = db_params$subcollection
      )
    },
    error = function(e) {
      msigdbr(
        species = species,
        category = db_params$collection,
        subcategory = db_params$subcollection
      )
    }
  )

  # Load msigdbr data
  msigdb_init <- msigdb_raw %>%
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

#' Build a pathway graph from custom gene sets and relations
#'
#' Convert user-supplied pathway gene sets and parent-child relations into the
#' internal objects used throughout the package.
#'
#' @details
#' Use this function when you have your own pathway collection and want to use
#' the plotting and Shiny interfaces without relying on the built-in Reactome or
#' GO resources.
#'
#' `gene_sets` must be a named list whose names are the human-readable pathway
#' names you want to use in plots. `pathways_relation` should contain pathway
#' identifiers used in the graph edges. If those identifiers differ from the
#' names in `gene_sets`, provide `translation_layer` to link them.
#'
#' @importFrom igraph graph_from_data_frame V vertex edges
#' @importFrom dplyr left_join
#' @importFrom stats setNames
#' @importFrom stringr str_to_title
#' @param gene_sets Named list of gene vectors. List names are treated as the
#'   pathway names shown to the user.
#' @param pathways_relation Data frame containing at least `parent` and `child`
#'   columns describing the directed pathway hierarchy.
#' @param translation_layer Optional two-column data frame that maps gene-set
#'   names to relation-layer identifiers. The first column should match names in
#'   `gene_sets`; the second should match the identifiers used in
#'   `pathways_relation`.
#' @return A list with `gene_sets`, `mapping`, and `pathway_graph` elements.
#' @examples
#' gene_sets <- list(
#'   Glycolysis = c("HK1", "PFKM", "PKM"),
#'   TCA = c("CS", "ACO2", "IDH3A")
#' )
#' pathways_relation <- data.frame(parent = "Glycolysis", child = "TCA")
#' prep <- load_and_preprocess_gene_sets(gene_sets, pathways_relation)
#' prep$mapping
#' @seealso [load_default_data()], [plot_subgraph()]
#' @export
load_and_preprocess_gene_sets <- function(gene_sets,
                                          pathways_relation,
                                          translation_layer = NULL) {
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a named list with pathway names as names.")
  }
  
  if (!is.data.frame(pathways_relation) || !all(c("parent", "child") %in% colnames(pathways_relation))) {
    stop("pathways_relation must be a data frame with 'parent' and 'child' columns.")
  }
  
  if (!is.null(translation_layer)) {
    if (!is.data.frame(translation_layer) || ncol(translation_layer) != 2) {
      stop("translation_layer must be a data frame with two columns.")
    }
    colnames(translation_layer) <- c("gene_sets_name", "relation_name")
  }
  
  gene_sets_names <- names(gene_sets)
  mapping <- data.frame(
    processed_name = gene_sets_names,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(translation_layer)) {
    mapping <- left_join(mapping,
                                translation_layer,
                                by = c("processed_name" = "gene_sets_name"))
    mapping$exact_source <- mapping$relation_name
  } else {
    mapping$exact_source <- mapping$processed_name
  }
  
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
    ifelse(grepl("REACTOME", complete_mapping$exact_source),
           str_to_title(gsub("_", " ", gsub("^REACTOME_", "", complete_mapping$exact_source))),
           str_to_title(gsub("_", " ", gsub("^GOBP_", "", complete_mapping$exact_source)))
    ),
    complete_mapping$processed_name
  )
  
  exact_source_to_name <- setNames(complete_mapping$processed_name, complete_mapping$exact_source)
  
  valid_exact_sources <- mapping$exact_source
  pathways_relation_filtered <- pathways_relation[pathways_relation$parent %in% valid_exact_sources &
                                                    pathways_relation$child %in% valid_exact_sources, ]
  
  g <- graph_from_data_frame(d = pathways_relation_filtered, directed = TRUE)
  
  V(g)$label <- exact_source_to_name[V(g)$name]
  missing_labels <- is.na(V(g)$label)
  V(g)$label[missing_labels] <- V(g)$name[missing_labels]
  
  V(g)$genes <- gene_sets[exact_source_to_name[V(g)$name]]
  V(g)$genes <- lapply(V(g)$genes,
                               function(x) if (is.null(x)) character(0) else x)
  
  return(list(
    gene_sets = gene_sets,
    mapping = mapping,
    pathway_graph = g
  ))
}


#' Add gene-set overlap information to graph edges
#'
#' Compute the number of shared genes between the source and target pathway of
#' every edge in a pathway graph and store the results as edge attributes.
#'
#' @details
#' The function adds two edge attributes to `g`: `overlap`, the number of shared
#' genes, and `percent_overlap`, the overlap as a percentage of the source
#' pathway gene-set size.
#'
#' @importFrom igraph E ends V
#' @param g Pathway graph as an `igraph` object. Vertices must already contain a
#'   `genes` attribute.
#' @param gene_sets Named list of gene sets. This is included for API symmetry
#'   with the rest of the package; overlap is derived from the vertex `genes`
#'   attributes on `g`.
#' @return The input graph with `overlap` and `percent_overlap` edge attributes
#'   added.
#' @examples
#' gene_sets <- list(Parent = c("G1", "G2"), Child = c("G2", "G3"))
#' prep <- load_and_preprocess_gene_sets(
#'   gene_sets,
#'   data.frame(parent = "Parent", child = "Child")
#' )
#' g2 <- add_overlap_to_edges(prep$pathway_graph, prep$gene_sets)
#' igraph::edge_attr(g2, "overlap")
#' @seealso [plot_subgraph()]
#' @export
add_overlap_to_edges <- function(g, gene_sets) {
  edge_ends <- ends(g, E(g))
  from_nodes <- edge_ends[, 1]
  to_nodes <- edge_ends[, 2]
  
  from_genes_list <- V(g)$genes[match(from_nodes, V(g)$name)]
  to_genes_list <- V(g)$genes[match(to_nodes, V(g)$name)]
  
  n_edges <- length(from_genes_list)
  overlaps <- vapply(seq_len(n_edges), function(i) {
    length(intersect(from_genes_list[[i]], to_genes_list[[i]]))
  }, FUN.VALUE = integer(1))
  
  from_sizes <- lengths(from_genes_list)
  percent_overlaps <- ifelse(from_sizes > 0, (overlaps / from_sizes) * 100, 0)
  
  E(g)$overlap <- overlaps
  E(g)$percent_overlap <- percent_overlaps
  
  return(g)
}

#' Get all descendant pathways below a parent pathway
#'
#' Return the processed pathway names for every descendant of a selected parent
#' pathway in the directed pathway graph.
#'
#' @importFrom igraph subcomponent V
#' @param parent_geneset_name Character scalar giving the parent pathway name in
#'   processed-name form.
#' @param mapping Mapping data frame containing `processed_name` and
#'   `exact_source` columns.
#' @param g Pathway graph as an `igraph` object.
#' @return Character vector of descendant pathway names, excluding the parent
#'   pathway itself.
#' @examples
#' gene_sets <- list(Parent = "G1", Child = "G2", Grandchild = "G3")
#' relations <- data.frame(
#'   parent = c("Parent", "Child"),
#'   child = c("Child", "Grandchild")
#' )
#' prep <- load_and_preprocess_gene_sets(gene_sets, relations)
#' get_child_pathways("Parent", prep$mapping, prep$pathway_graph)
#' @seealso [plot_subgraph()], [load_and_preprocess_gene_sets()]
#' @export
get_child_pathways <- function(parent_geneset_name, mapping, g) {
  if (length(parent_geneset_name) == 0) {
    stop("No parent_geneset_name provided or it is an empty vector.")
  }
  
  parent_indices <- match(parent_geneset_name, mapping$processed_name)
  
  if (length(parent_indices) == 0 || is.na(parent_indices)) {
    stop("Parent pathway not found in the mapping.")
  }
  
  parent_exact_source <- mapping$exact_source[parent_indices]
  
  if (length(parent_exact_source) == 0 || is.na(parent_exact_source)) {
    stop("Parent pathway not found in the mapping.")
  }
  
  parent_vertex <- which(V(g)$name == parent_exact_source)
  if (length(parent_vertex) == 0) {
    stop("Parent pathway not found in the graph.")
  }
  
  descendants <- subcomponent(g, v = parent_vertex, mode = "out")
  
  descendant_exact_sources <- V(g)$name[descendants]
  descendant_exact_sources <- descendant_exact_sources[descendant_exact_sources != parent_exact_source]
  
  descendant_geneset_names <- mapping$processed_name[match(descendant_exact_sources, mapping$exact_source)]
  descendant_geneset_names <- descendant_geneset_names[!is.na(descendant_geneset_names)]
  
  return(descendant_geneset_names)
}
