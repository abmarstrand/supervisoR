#' Plot a pathway subgraph
#'
#' Visualize a pathway and all of its descendants as a graph. When multiple
#' conditions are supplied, nodes are rendered as pathway glyphs summarizing the
#' enrichment profile across conditions. When exactly one condition is supplied,
#' nodes are drawn as colored circles whose size scales with the absolute score.
#'
#' @details
#' This is the main static plotting function in the package. It can work either
#' with custom graph objects supplied through `mapping`, `g`, and `gene_sets`, or
#' it can lazily load the built-in Reactome/GO resources through
#' [load_default_data()].
#'
#' `enrichment_scores` may be provided in one of two ways:
#'
#' * pathway names already stored in `rownames(enrichment_scores)` and
#'   `pathway_column = NULL`
#' * pathway names stored in one column of `enrichment_scores`, in which case
#'   `pathway_column` must point to that column
#'
#' Pathway names are matched against `mapping$processed_name`, not the internal
#' `exact_source` identifiers.
#'
#' @importFrom igraph subcomponent induced_subgraph V delete_vertices distances edge_attr neighbors add_vertices add_edges E V<- E<-
#' @importFrom ggraph ggraph create_layout geom_edge_link scale_edge_color_manual scale_edge_width geom_node_point circle scale_edge_linetype_manual geom_node_text
#' @importFrom ggplot2 theme_void ggplotGrob geom_label ylab ggtitle
#' @importFrom cowplot draw_grob ggdraw plot_grid
#' @importFrom ggimage geom_image
#' @importFrom ggrepel geom_label_repel
#' @importFrom colorspace scale_colour_continuous_diverging
#' @importFrom stringr str_wrap
#' @importFrom grid arrow unit
#' @param parent_geneset_name Character scalar giving the pathway to plot. This
#'   should use the processed pathway names shown in the mapping table.
#' @param enrichment_scores Optional data frame of enrichment scores. Rows must
#'   correspond to pathways and columns to conditions after any `pathway_column`
#'   processing.
#' @param conditions Optional character vector of conditions to plot. If `NULL`,
#'   all available condition columns are used.
#' @param significance_scores Optional data frame of significance values aligned
#'   with `enrichment_scores`.
#' @param pathway_column Column index, column name, or `NULL` describing where
#'   pathway names are stored in `enrichment_scores` and `significance_scores`.
#'   Use `NULL` when row names already contain pathway names.
#' @param species Species passed to [load_default_data()] when default resources
#'   are loaded automatically.
#' @param database Database passed to [load_default_data()]. Supported built-in
#'   values are currently `"reactome"` and `"GO"`.
#' @param mapping Optional mapping data frame with `processed_name` and
#'   `exact_source` columns.
#' @param g Optional pathway graph as an `igraph` object.
#' @param gene_sets Optional named list of pathway gene sets.
#' @param layout Layout algorithm passed to `ggraph::create_layout()`.
#' @param circular Logical; if `TRUE`, request a circular layout.
#' @param hide_nodes_without_enrichment Logical; if `TRUE`, remove nodes for
#'   which all selected enrichment values are missing.
#' @param enrichment_limits Optional numeric vector of length 2 controlling the
#'   color scale and y-axis limits. When `NULL`, limits are inferred from the
#'   visible scores.
#' @param use_node_label Logical; if `TRUE`, pathway labels are drawn next to the
#'   nodes.
#' @param reference Optional condition name to subtract from every selected
#'   condition before plotting.
#' @param adjust_edge_thickness Logical; if `TRUE`, scale edge thickness by gene
#'   set overlap.
#' @param edge_percentage_labels Logical; if `TRUE`, annotate edges with the
#'   percentage overlap relative to the parent/source pathway gene set.
#' @param glyph_size Numeric vector of length 2 giving glyph width and height in
#'   pixels.
#' @param res Numeric rendering resolution in dpi.
#' @param legend_position Legend position passed into the assembled plot. Use
#'   `"none"` to suppress legends.
#' @param label_wrap_width Integer width used to wrap node labels.
#' @param edge_arrows Logical; if `TRUE`, draw arrowheads on edges.
#' @param repel_labels Logical; if `TRUE`, use `ggrepel` for node labels.
#' @param node_label_size Numeric label size.
#' @param max_depth Optional non-negative integer. When supplied, the plot is
#'   truncated to descendants within that graph distance and ellipsis nodes are
#'   added where deeper descendants were removed.
#' @param lollipop_plot Logical; if `TRUE`, multi-condition glyphs use a lollipop
#'   representation. Otherwise they use a bar representation.
#' @param lollipop_colors Optional named character vector mapping condition names
#'   to colors.
#' @param lollipop_significance Logical; if `TRUE`, encode significance inside
#'   lollipop glyphs.
#' @param lollipop_size_scale Numeric scalar that scales both the lollipop stem
#'   width and point size inside pathway glyphs. Values above `1` enlarge the
#'   lollipops and values below `1` shrink them.
#' @param significance_cutoff Numeric significance threshold used when
#'   `lollipop_significance = TRUE`.
#' @param repel Deprecated alias for `repel_labels`.
#' @return A `ggplot` object representing the selected pathway subgraph.
#' @examples
#' gene_sets <- list(
#'   Parent = c("G1", "G2", "G3"),
#'   ChildA = c("G1", "G2"),
#'   ChildB = c("G3", "G4")
#' )
#' relations <- data.frame(
#'   parent = c("Parent", "Parent"),
#'   child = c("ChildA", "ChildB")
#' )
#' prep <- load_and_preprocess_gene_sets(gene_sets, relations)
#' scores <- data.frame(
#'   cond1 = c(1.2, 0.4, -0.8),
#'   cond2 = c(0.7, -0.2, -1.1),
#'   row.names = c("Parent", "ChildA", "ChildB")
#' )
#' plot_subgraph(
#'   parent_geneset_name = "Parent",
#'   enrichment_scores = scores,
#'   conditions = c("cond1", "cond2"),
#'   pathway_column = NULL,
#'   mapping = prep$mapping,
#'   g = prep$pathway_graph,
#'   gene_sets = prep$gene_sets,
#'   lollipop_plot = TRUE
#' )
#' @seealso [load_default_data()], [load_and_preprocess_gene_sets()],
#'   [run_pathway_shiny_app()]
#' @export
plot_subgraph <- function(parent_geneset_name, enrichment_scores = NULL, conditions = NULL,
                          significance_scores = NULL,
                          pathway_column = 1L, species = "Homo sapiens", database = "reactome",
                          mapping = NULL, g = NULL, gene_sets = NULL,
                          layout = "kk", circular = FALSE,
                          hide_nodes_without_enrichment = FALSE,
                          enrichment_limits = NULL, use_node_label = TRUE,
                          reference = NULL, adjust_edge_thickness = FALSE,
                          edge_percentage_labels = FALSE, glyph_size = c(80, 60),
                          res = 96, legend_position = NULL,
                          label_wrap_width = 40, edge_arrows = TRUE,
                          repel_labels = FALSE, node_label_size = 1, max_depth = NULL,
                          lollipop_plot = FALSE, lollipop_colors = NULL,
                          lollipop_significance = FALSE,
                          significance_cutoff = 0.05,
                          lollipop_size_scale = 1,
                          repel = NULL) {
  if (!is.null(repel)) {
    repel_labels <- isTRUE(repel)
  }
  
  if (is.null(mapping) || is.null(g) || is.null(gene_sets)) {
    data <- load_default_data(database = database, species = species)
    gs <- data$gene_sets
    rel <- data$relations
    trans <- data.frame(
      gene_sets_name = data$mapping$processed_name,
      relation_name = data$mapping$exact_source,
      stringsAsFactors = FALSE
    )
    data <- load_and_preprocess_gene_sets(
      gene_sets = gs,
      pathways_relation = rel,
      translation_layer = trans
    )
    if (is.null(mapping)) mapping <- data$mapping
    if (is.null(g)) g <- data$pathway_graph
    if (is.null(gene_sets)) gene_sets <- data$gene_sets
  }
  
  if (!is.null(enrichment_scores) && !is.null(pathway_column)) {
    if (is.character(enrichment_scores[, pathway_column])) {
      rownames(enrichment_scores) <- enrichment_scores[, pathway_column]
      enrichment_scores <- enrichment_scores[, -pathway_column, drop = FALSE]
    } else {
      stop("Pathway column does not contain characters. Check whether you selected the correct pathway column.")
    }
  }
  
  if (!is.null(significance_scores) && !is.null(pathway_column)) {
    if (is.character(significance_scores[, pathway_column])) {
      rownames(significance_scores) <- significance_scores[, pathway_column]
      significance_scores <- significance_scores[, -pathway_column, drop = FALSE]
    } else {
      stop("Pathway column in significance_scores does not contain characters.")
    }
  }
  
  if (!parent_geneset_name %in% mapping$processed_name) {
    stop("Parent pathway not found in the mapping.")
  }
  
  parent_exact_source <- mapping$exact_source[match(parent_geneset_name, mapping$processed_name)]
  parent_vertex <- which(V(g)$name == parent_exact_source)
  if (length(parent_vertex) == 0) {
    stop("Parent pathway not found in the graph.")
  }
  
  if (is.null(max_depth)) {
    descendants <- subcomponent(g, v = parent_vertex, mode = "out")
    sub_g <- induced_subgraph(g, vids = descendants)
  } else {
    dist_vec <- distances(g, v = parent_vertex, mode = "out")[1, ]
    allowed_nodes <- names(dist_vec)[dist_vec <= max_depth]
    sub_g <- induced_subgraph(g, vids = allowed_nodes)
    border_nodes <- names(dist_vec)[dist_vec == max_depth]
    for (bn in border_nodes) {
      bn_idx <- which(V(g)$name == bn)
      out_neighbors <- neighbors(g, v = bn_idx, mode = "out")
      out_neighbor_names <- V(g)$name[out_neighbors]
      truncated <- out_neighbor_names[!(out_neighbor_names %in% allowed_nodes)]
      if (length(truncated) > 0) {
        ellipsis_node_name <- paste0("ellipsis_for_", bn)
        while (ellipsis_node_name %in% V(sub_g)$name) {
          ellipsis_node_name <- paste0(ellipsis_node_name, "_x")
        }
        sub_g <- add_vertices(sub_g, nv = 1, name = ellipsis_node_name)
        sub_g <- add_edges(sub_g, c(bn, ellipsis_node_name), relation = "truncated")
        V(sub_g)$label[V(sub_g)$name == ellipsis_node_name] <- "..."
      }
    }
  }
  
  normal_nodes <- V(sub_g)$name[!grepl("^ellipsis_for_", V(sub_g)$name)]
  V(sub_g)$label[!grepl("^ellipsis_for_", V(sub_g)$name)] <-
    mapping$processed_name[match(normal_nodes, mapping$exact_source)]
  
  missing_labels <- is.na(V(sub_g)$label) &
    !grepl("^ellipsis_for_", V(sub_g)$name) &
    V(sub_g)$name != parent_exact_source
  if (any(missing_labels)) {
    sub_g <- delete_vertices(sub_g, V(sub_g)[missing_labels])
  }
  
  if (!is.null(enrichment_scores)) {
    conditions <- resolve_conditions(enrichment_scores, conditions)
    if (!is.null(reference) && reference %in% colnames(enrichment_scores)) {
      enrichment_scores <- sweep(enrichment_scores, 1, enrichment_scores[, reference], FUN = "-")
      conditions <- setdiff(conditions, reference)
    }
    if (!is.null(significance_scores)) {
      conditions <- resolve_conditions(significance_scores, conditions)
    }
  }
  
  if (!is.null(enrichment_scores) && length(conditions) == 0) {
    stop("No conditions remain after applying the reference condition.")
  }
  
  if (!is.null(enrichment_scores) && !is.null(conditions)) {
    node_labels <- V(sub_g)$label
    valid_labels <- node_labels[node_labels %in% rownames(enrichment_scores)]
    scores <- if (length(valid_labels) > 0) {
      enrichment_scores[valid_labels, conditions, drop = FALSE]
    } else {
      NULL
    }
    
    missing_scores <- if (length(node_labels) > 0) {
      tmp <- enrichment_scores[node_labels, conditions, drop = FALSE]
      rowSums(is.na(tmp)) == length(conditions)
    } else {
      logical(0)
    }
    missing_scores[is.na(missing_scores)] <- TRUE
    
    if (hide_nodes_without_enrichment && any(missing_scores)) {
      keep_nodes <- !missing_scores | grepl("^ellipsis_for_", V(sub_g)$name)
      sub_g <- induced_subgraph(sub_g, vids = V(sub_g)[keep_nodes])
      valid_labels <- V(sub_g)$label[V(sub_g)$label %in% rownames(enrichment_scores)]
      scores <- if (length(valid_labels) > 0) {
        enrichment_scores[valid_labels, conditions, drop = FALSE]
      } else {
        NULL
      }
    }
    
    enrichment_limits <- normalize_enrichment_limits(
      enrichment_scores = scores,
      enrichment_limits = enrichment_limits,
      zero_anchor = TRUE,
      symmetric = TRUE,
      padding = 0.08
    )
  }
  
  if (adjust_edge_thickness || edge_percentage_labels) {
    sub_g <- add_overlap_to_edges(sub_g, gene_sets)
  }
  
  layout_data <- create_layout(sub_g, layout = layout, circular = circular)
  layout_data$wrapped_label <- str_wrap(layout_data$label, width = label_wrap_width)
  is_ellipsis <- layout_data$label == "..."
  single_condition <- !is.null(conditions) && length(conditions) == 1
  
  if (!single_condition && !is.null(conditions) && !is.null(enrichment_scores)) {
    all_labels <- V(sub_g)$label
    glyphs <- lapply(all_labels, function(pathway) {
      if (is.na(pathway) || identical(pathway, "...") || !(pathway %in% rownames(enrichment_scores))) {
        return(NA_character_)
      }
      pathway_scores <- enrichment_scores[pathway, conditions, drop = FALSE]
      if (all(is.na(pathway_scores))) {
        return(NA_character_)
      }
      pathway_significance <- NULL
      if (!is.null(significance_scores) && pathway %in% rownames(significance_scores)) {
        pathway_significance <- significance_scores[pathway, conditions, drop = FALSE]
      }
      img <- create_glyph_on_the_fly(
        pathway = pathway,
        conditions = conditions,
        enrichment_scores = enrichment_scores,
        enrichment_limits = enrichment_limits,
        significance_scores = pathway_significance,
        glyph_size = glyph_size,
        res = res,
        lollipop_plot = lollipop_plot,
        lollipop_colors = lollipop_colors,
        lollipop_significance = lollipop_significance,
        significance_cutoff = significance_cutoff,
        lollipop_size_scale = lollipop_size_scale
      )
      if (is.null(img)) NA_character_ else img[[2]]
    })
    layout_data$image <- vapply(glyphs, identity, FUN.VALUE = character(1))
  } else if (single_condition && !is.null(conditions) && !is.null(enrichment_scores)) {
    layout_data$enrichment_score <- enrichment_scores[layout_data$label, conditions[1]]
    abs_score <- abs(layout_data$enrichment_score)
    abs_score[!is.finite(abs_score)] <- 0
    score_span <- max(abs(enrichment_limits), na.rm = TRUE)
    if (!is.finite(score_span) || score_span <= 0) {
      score_span <- 1
    }
    layout_data$node_size <- 3 + 5 * (abs_score / score_span)
  }
  
  relation_linetypes <- c(
    "is_a" = "solid",
    "part_of" = "dashed",
    "regulates" = "dotted",
    "positively_regulates" = "dotdash",
    "negatively_regulates" = "longdash",
    "truncated" = "solid"
  )
  relation_colors <- c(
    "is_a" = "black",
    "part_of" = "darkcyan",
    "regulates" = "orange",
    "positively_regulates" = "steelblue",
    "negatively_regulates" = "salmon",
    "truncated" = "grey65"
  )
  
  E(sub_g)$relation <- factor(E(sub_g)$relation, levels = names(relation_linetypes))
  
  p <- ggraph(layout_data) +
    geom_edge_link(
      aes(
        width = if (adjust_edge_thickness) percent_overlap else 1,
        label = if (edge_percentage_labels) paste0(round(percent_overlap, 1), "%") else NA,
        color = relation,
        linetype = relation
      ),
      show.legend = FALSE,
      arrow = if (edge_arrows) arrow(length = unit(0.0075, "npc"), type = "open") else NULL,
      end_cap = if (edge_arrows) circle(0.015, "npc") else NULL,
      start_cap = if (edge_arrows) circle(0.015, "npc") else NULL,
      angle_calc = "along",
      check_overlap = TRUE,
      label_dodge = unit(0.015, "npc"),
      label_push = unit(0.015, "npc"),
      alpha = 0.7
    ) +
    scale_edge_width(range = c(0.4, 2.4)) +
    scale_edge_color_manual(values = relation_colors, na.value = "black") +
    scale_edge_linetype_manual(values = relation_linetypes, na.value = "solid") +
    theme_void()
  
  has_image <- if ("image" %in% names(layout_data)) !is.na(layout_data$image) else rep(FALSE, nrow(layout_data))
  
  if (!single_condition && !is.null(conditions) && !is.null(enrichment_scores)) {
    num_nodes <- max(nrow(layout_data), 1L)
    glyph_size_scale <- min(0.11, max(0.03, 0.55 / sqrt(num_nodes)))
    if (any(has_image & !is_ellipsis)) {
      p <- p + geom_image(
        data = subset(layout_data, has_image & !is_ellipsis),
        aes(x = x, y = y, image = image),
        size = glyph_size_scale,
        by = "width",
        asp = 1.4,
        na.rm = TRUE
      )
    }
    if (any(!has_image & !is_ellipsis)) {
      p <- p + geom_node_point(
        data = subset(layout_data, !has_image & !is_ellipsis),
        aes(x = x, y = y),
        color = "#6A3D9A",
        size = glyph_size_scale * 22,
        alpha = 0.9
      )
    }
  } else if (single_condition && !is.null(conditions) && !is.null(enrichment_scores)) {
    p <- p +
      geom_node_point(
        data = subset(layout_data, !is_ellipsis),
        aes(x = x, y = y, color = enrichment_score, size = node_size),
        alpha = 0.95
      ) +
      scale_colour_continuous_diverging(palette = "Berlin", limits = enrichment_limits) +
      ggplot2::scale_size_continuous(range = c(3, 8), guide = "none")
  } else {
    p <- p + geom_node_point(
      data = subset(layout_data, !is_ellipsis),
      aes(x = x, y = y),
      size = 4.5,
      color = "steelblue"
    )
  }
  
  if (any(is_ellipsis)) {
    p <- p + geom_node_text(
      data = subset(layout_data, is_ellipsis),
      aes(x = x, y = y, label = label),
      size = 4,
      color = "grey20",
      fontface = "bold"
    )
  }
  
  if (use_node_label && any(!is_ellipsis)) {
    non_ellipsis_data <- subset(layout_data, !is_ellipsis)
    if (repel_labels) {
      p <- p + geom_label_repel(
        data = non_ellipsis_data,
        aes(x = x, y = y, label = wrapped_label),
        size = 3 * node_label_size,
        label.size = 0.15,
        label.padding = unit(0.12, "lines"),
        fill = "white",
        box.padding = unit(0.4, "lines"),
        point.padding = unit(0.35, "lines"),
        force = 2
      )
    } else {
      y_range <- diff(range(non_ellipsis_data$y, na.rm = TRUE))
      if (!is.finite(y_range) || y_range == 0) {
        y_range <- 1
      }
      gs_val <- if (exists("glyph_size_scale")) glyph_size_scale else 0.05
      nudge_y <- -max(gs_val * y_range * 0.35, 0.08)
      p <- p + geom_label(
        data = non_ellipsis_data,
        aes(x = x, y = y, label = wrapped_label),
        size = 3 * node_label_size,
        vjust = 1,
        nudge_y = nudge_y,
        label.padding = unit(0.1, "lines"),
        label.size = 0.15,
        fill = "white"
      )
    }
  }
  
  node_x <- layout_data$x
  node_y <- layout_data$y
  corners <- data.frame(
    corner = c("topright", "topleft", "bottomleft", "bottomright", "none"),
    x = c(max(node_x), min(node_x), min(node_x), max(node_x), max(node_x)),
    y = c(max(node_y), max(node_y), min(node_y), min(node_y), max(node_y))
  )
  if (!is.null(legend_position) && legend_position %in% corners$corner) {
    best_corner <- corners[corners$corner == legend_position, , drop = FALSE]
  } else {
    corners$min_dist <- vapply(seq_len(nrow(corners)), function(i) {
      sqrt(min((node_x - corners$x[i])^2 + (node_y - corners$y[i])^2))
    }, FUN.VALUE = numeric(1))
    best_corner <- corners[which.max(corners$min_dist), , drop = FALSE]
  }
  
  if (!single_condition && !is.null(enrichment_scores) && !is.null(conditions)) {
    legend_plot <- generate_legend_plot(
      conditions = conditions,
      enrichment_limits = enrichment_limits,
      reference = reference,
      lollipop_plot = lollipop_plot,
      lollipop_colors = lollipop_colors,
      lollipop_significance = lollipop_significance,
      lollipop_size_scale = lollipop_size_scale
    )
    relations <- if ("relation" %in% names(edge_attr(sub_g))) unique(E(sub_g)$relation) else character(0)
    relations_present <- relations[!is.na(relations)]
    relations_present <- relations_present[relations_present %in% names(relation_colors)]
    if (length(relations_present) > 0) {
      relation_legend <- generate_relation_legend(
        relations_present = relations_present,
        relation_colors = relation_colors,
        relation_linetypes = relation_linetypes,
        textsize = 3
      )
      combined_plot <- plot_grid(legend_plot, relation_legend, ncol = 1, rel_heights = c(3, 1))
      legend_grob <- ggplotGrob(combined_plot)
    } else {
      legend_grob <- ggplotGrob(legend_plot)
    }
    
    legend_coords <- switch(
      best_corner$corner,
      "topleft" = list(x = 0.03, y = 0.73),
      "topright" = list(x = 0.74, y = 0.73),
      "bottomleft" = list(x = 0.03, y = 0.03),
      "bottomright" = list(x = 0.74, y = 0.03),
      "none" = list(x = NA, y = NA),
      list(x = 0.74, y = 0.03)
    )
    if (!identical(best_corner$corner, "none")) {
      p <- ggdraw(p) +
        draw_grob(legend_grob, x = legend_coords$x, y = legend_coords$y, width = 0.23, height = 0.24)
    }
  }
  
  p
}
