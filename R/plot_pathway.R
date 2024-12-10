#' Plot a Subgraph with Glyphs as Nodes
#'
#' Plots a subgraph with nodes represented by glyphs that visualize enrichment scores across specified conditions.
#' Glyphs are generated on-the-fly using the \code{create_glyph_on_the_fly} function, allowing dynamic comparisons
#' and efficient plotting even with a large number of conditions.
#'
#' @importFrom igraph subcomponent induced_subgraph V delete_vertices ends `V<-` E `E<-` distances edge_attr neighbors
#' @importFrom ggraph ggraph create_layout geom_edge_link geom_node_label geom_node_text scale_edge_color_manual  scale_edge_width geom_node_point circle scale_edge_linetype_manual
#' @importFrom ggplot2 theme_void ggplotGrob geom_label expansion ylim xlim xlab ylab ggtitle
#' @importFrom cowplot draw_grob ggdraw plot_grid
#' @importFrom ggimage geom_image
#' @importFrom ggrepel geom_label_repel
#' @importFrom magrittr %>%
#' @importFrom colorspace scale_colour_continuous_diverging
#' @importFrom stringr str_wrap
#' @importFrom grid arrow unit
#' @param parent_geneset_name Name of the parent geneset. This is used to extract the subgraph of interest.
#' @param enrichment_scores Data frame of enrichment scores (rows are pathways, columns are conditions).
#' @param conditions Character vector of conditions to include in the glyphs.
#' @param mapping Data frame containing mapping information.
#' @param g The pathway graph (igraph object).
#' @param gene_sets List of gene sets used for edge calculations.
#' @param layout The layout algorithm to use (e.g., "kk", "fr"). Default is "kk".
#' @param circular Logical, whether to use a circular layout. Default is \code{FALSE}.
#' @param hide_nodes_without_enrichment Logical, whether to hide nodes without enrichment scores. Default is \code{TRUE}.
#' @param enrichment_limits Numeric vector of length 2 specifying the enrichment limits. If \code{NULL}, limits are calculated from the data.
#' @param use_node_label Logical, whether to display node labels. Default is \code{TRUE}.
#' @param reference Optional character string specifying the reference condition.
#' @param adjust_edge_thickness Logical, whether to adjust edge thickness based on overlap. Default is \code{FALSE}.
#' @param edge_percentage_labels Logical, whether to add edge percentage labels. Default is \code{FALSE}.
#' @param glyph_size Numeric vector specifying the width and height of the glyphs in pixels. Default is \code{c(80, 60)}.
#' @param res Numeric value specifying the resolution of the glyph images in dpi. Default is 96.
#' @param legend_position Optional character string ("topright", "topleft", "bottomleft" or "bottomright") determining final legend position.
#' @param label_wrap_width Integer specifying the maximum number of characters in each line of the label. Default is 40.
#' @param edge_arrows Logical, whether to display edges as arrows pointing from parents to children. Default is \code{TRUE}.
#' @param repel_labels Logical, whether to repel labels to avoid overlaps. Default is \code{FALSE}. If you have a very large graph it is reccommended to try this out.
#' @param max_depth Integer, specifying the maximal depth shown at one time on the graph. Defaults to \code{NULL}. If you have a very large graph it is reccommended to try this out.
#' @return A ggplot object representing the subgraph.
#' @export
plot_subgraph <- function(parent_geneset_name, enrichment_scores = NULL, conditions = NULL,
                          mapping, g, gene_sets = NULL,
                          layout = "kk", circular = FALSE,
                          hide_nodes_without_enrichment = TRUE,
                          enrichment_limits = NULL, use_node_label = TRUE,
                          reference = NULL, adjust_edge_thickness = FALSE,
                          edge_percentage_labels = FALSE, glyph_size = c(80, 60),
                          res = 96, legend_position = NULL,
                          label_wrap_width = 40, edge_arrows = TRUE,
                          repel_labels = FALSE, max_depth = NULL) {
  if (!parent_geneset_name %in% mapping$processed_name) {
    stop("Parent pathway not found in the mapping.")
  }

  parent_exact_source <- mapping$exact_source[match(parent_geneset_name, mapping$processed_name)]

  parent_vertex <- which(V(g)$name == parent_exact_source)
  if (length(parent_vertex) == 0) {
    stop("Parent pathway not found in the graph.")
  }

  # If max_depth is NULL, get all descendants as before
  if (is.null(max_depth)) {
    descendants <- subcomponent(g, v = parent_vertex, mode = "out")
    sub_g <- induced_subgraph(g, vids = descendants)
  } else {
    # Compute distances
    dist_vec <- distances(g, v = parent_vertex, mode = "out")[1, ]
    allowed_nodes <- names(dist_vec)[dist_vec <= max_depth]
    sub_g <- induced_subgraph(g, vids = allowed_nodes)
    # Identify border nodes
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

  normal_nodes <- V(sub_g)$name[!(grepl("^ellipsis_for_", V(sub_g)$name))]
  V(sub_g)$label[!(grepl("^ellipsis_for_", V(sub_g)$name))] <- mapping$processed_name[match(normal_nodes, mapping$exact_source)]

  missing_labels <- is.na(V(sub_g)$label) & !grepl("^ellipsis_for_", V(sub_g)$name) & V(sub_g)$name != parent_exact_source
  if (any(missing_labels)) {
    sub_g <- delete_vertices(sub_g, V(sub_g)[missing_labels])
  }

  if (!is.null(enrichment_scores) && !is.null(conditions)) {
    if (!is.null(reference) && reference %in% conditions) {
      enrichment_scores <- enrichment_scores - enrichment_scores[, reference, drop = TRUE]
      conditions <- setdiff(conditions, reference)
    }
    node_labels <- V(sub_g)$label
    valid_labels <- node_labels[node_labels %in% rownames(enrichment_scores)]
    scores <- enrichment_scores[valid_labels, conditions, drop = FALSE]
    missing_scores <- rowSums(is.na(enrichment_scores[node_labels,
                                                      conditions,
                                                      drop = FALSE])) == length(conditions)
    missing_scores[is.na(missing_scores)] <- TRUE
    if (hide_nodes_without_enrichment && any(missing_scores)) {
      keep_nodes <- !missing_scores | grepl("^ellipsis_for_", V(sub_g)$name)
      sub_g <- induced_subgraph(sub_g, vids = V(sub_g)[keep_nodes])
      valid_labels <- V(sub_g)$label[V(sub_g)$label %in% valid_labels]
      if (length(valid_labels) > 0) {
        scores <- enrichment_scores[valid_labels, conditions, drop = FALSE]
      } else {
        scores <- NULL
      }
    }
    if (is.null(enrichment_limits) && !is.null(scores)) {
      enrichment_limits <- range(scores, na.rm = TRUE)
    }
  }

  if (adjust_edge_thickness || edge_percentage_labels) {
    sub_g <- add_overlap_to_edges(sub_g, gene_sets)
  }

  layout_data <- create_layout(sub_g, layout = layout, circular = circular)
  layout_data$wrapped_label <- str_wrap(layout_data$label, width = label_wrap_width)
  single_condition <- !is.null(conditions) && length(conditions) == 1

  if (!single_condition && !is.null(conditions) && !is.null(enrichment_scores)) {
    glyphs <- lapply(V(sub_g)$label, function(pathway) {
      if (pathway == "...") return(NA)
      pathway_scores <- enrichment_scores[pathway, conditions, drop = FALSE]
      img <- create_glyph_on_the_fly(
        pathway = pathway,
        conditions = conditions,
        enrichment_scores = pathway_scores,
        enrichment_limits = enrichment_limits,
        glyph_size = glyph_size,
        res = res
      )
      if (!is.null(img)) img[[2]] else NA
    })
    node_images <- sapply(glyphs, function(x) x)
    layout_data$image <- node_images
  } else if (single_condition && !is.null(conditions) && !is.null(enrichment_scores)) {
    single_condition_name <- conditions[1]
    layout_data$enrichment_score <- enrichment_scores[layout_data$label, single_condition_name]
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
    "truncated" = "grey"
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
      arrow = if (edge_arrows) grid::arrow(length = unit(0.0075, "npc"), type = "open") else NULL,
      end_cap = if (edge_arrows) circle(0.015, "npc") else NULL,
      start_cap = if (edge_arrows) circle(0.015, "npc") else NULL,
      angle_calc = "along",
      check_overlap = TRUE,
      label_dodge = unit(0.015, "npc"),
      label_push = unit(0.015, "npc")
    ) +
    scale_edge_width(range = c(0.5, 2)) +
    theme_void() +
    scale_edge_color_manual(values = relation_colors, na.value = "black") +
    scale_edge_linetype_manual(values = relation_linetypes, na.value = "solid")

  # Drawing nodes:
  # Multiple conditions with images
  if (!single_condition && !is.null(conditions) && !is.null(enrichment_scores)) {
    num_nodes <- nrow(layout_data)
    min_glyph_size <- 0.02
    max_glyph_size <- 0.1
    glyph_size_scale <- min(max_glyph_size, max(min_glyph_size, 0.5 / sqrt(num_nodes)))
    if (any(!is.na(layout_data$image) & layout_data$label != "...")) {
      p <- p + geom_image(data = subset(layout_data, !is.na(image) & label != "..."),
                          aes(x = x, y = y, image = image),
                          size = glyph_size_scale, na.rm = TRUE)
    }
    if (any(is.na(layout_data$image) & layout_data$label != "...")) {
      p <- p + geom_node_point(data = subset(layout_data, is.na(image) & label != "..."),
                               aes(x = x, y = y),
                               color = "purple", size = glyph_size_scale * 25)
    }
    # Ellipsis nodes: just text
    if (any(layout_data$label == "...")) {
      p <- p + geom_node_text(data = subset(layout_data, label == "..."),
                              aes(x = x, y = y, label = label),
                              size = 10, color = "black", fontface = "bold")
    }
  } else if (single_condition && !is.null(conditions) && !is.null(enrichment_scores)) {
    # Exclude ellipses from points
    p <- p + geom_node_point(data = subset(layout_data, label != "..."),
                             aes(color = enrichment_score), size = 5) +
      scale_colour_continuous_diverging(palette = "Berlin", limits = enrichment_limits)
    # Ellipsis nodes
    if (any(layout_data$label == "...")) {
      p <- p + geom_node_text(data = subset(layout_data, label == "..."),
                              aes(x = x, y = y, label = label),
                              size = 10, color = "black", fontface = "bold")
    }
  } else {
    # No enrichment scenario
    # Exclude ellipses from points
    p <- p + geom_node_point(data = subset(layout_data, label != "..."),
                             aes(x = x, y = y),
                             size = 5, color = "steelblue")
    # Ellipsis nodes
    if (any(layout_data$label == "...")) {
      p <- p + geom_node_text(data = subset(layout_data, label == "..."),
                              aes(x = x, y = y, label = label),
                              size = 10, color = "black", fontface = "bold")
    }
  }
  # Labels for non-ellipsis nodes if requested
  if (use_node_label && any(layout_data$label != "...")) {
    non_ellipsis_data <- subset(layout_data, label != "...")
    if (repel_labels) {
      p <- p + geom_label_repel(
        data = non_ellipsis_data,
        aes(x = x, y = y, label = wrapped_label),
        size = 3,
        box.padding = unit(0.5, "lines"),
        point.padding = unit(0.5, "lines"),
        force = 2
      )
    } else {
      if ("image" %in% names(layout_data)) {
        y_range <- diff(range(non_ellipsis_data$y))
        glyph_size_scale <- ifelse(exists("glyph_size_scale"), glyph_size_scale, 0.05)
        nudge_y <- -glyph_size_scale * y_range * 0.35
        p <- p + geom_label(
          data = non_ellipsis_data,
          aes(x = x, y = y, label = wrapped_label),
          size = 3,
          vjust = 1,
          nudge_y = nudge_y,
          label.padding = unit(0.1, "lines"),
          label.size = 0.2
        )
      } else {
        nudge_y <- -0.1
        p <- p + geom_label(
          data = non_ellipsis_data,
          aes(x = x, y = y, label = wrapped_label),
          size = 3,
          vjust = 1,
          nudge_y = nudge_y,
          label.padding = unit(0.1, "lines"),
          label.size = 0.2
        )
      }
    }
  }
  # Add legend if needed
  # Extract node positions
  node_x <- layout_data$x
  node_y <- layout_data$y
  # Determine plot boundaries
  x_min <- min(node_x)
  x_max <- max(node_x)
  y_min <- min(node_y)
  y_max <- max(node_y)

  corners <- data.frame(
    corner = c("topright", "topleft", "bottomleft", "bottomright"),
    x = c(x_max, x_min, x_min, x_max),
    y = c(y_max, y_max, y_min, y_min)
  )

  if (!is.null(legend_position) && legend_position %in% corners$corner) {
    best_corner <- corners[corners$corner == legend_position, ]
  } else {
    # Calculate minimum distance from each corner to any node
    corners$min_dist <- sapply(seq_len(nrow(corners)), function(i) {
      sqrt(min((node_x - corners$x[i])^2 + (node_y - corners$y[i])^2))
    })
    # Choose the corner with the maximum minimum distance
    best_corner <- corners[which.max(corners$min_dist), ]
  }

  # Create the appropriate legend based on the condition
  if (!single_condition && !is.null(enrichment_scores)) {
    legend_plot <- generate_legend_plot(
      conditions = conditions,
      enrichment_limits = enrichment_limits,
      reference = reference
    )
    if ("relation" %in% names(edge_attr(sub_g))) {
      relations <- unique(E(sub_g)$relation)
    } else {
      relations <- character(0)
    }
    relations_present <- relations[!is.na(relations)]
    relations_present <- relations_present[relations_present %in% names(relation_colors)]
    if (length(relations_present) > 0) {
      relation_legend <- generate_relation_legend(
        relations_present = relations_present,
        relation_colors = relation_colors,
        relation_linetypes = relation_linetypes,
        textsize = 3
      )
      # Combine the two legends vertically
      combined_plot <- plot_grid(legend_plot, relation_legend, ncol = 1, rel_heights = c(3, 1))
      legend_grob <- ggplotGrob(combined_plot)
    } else {
      legend_grob <- ggplotGrob(legend_plot)
    }
    legend_coords <- switch(best_corner$corner,
                            "topleft" = list(x = 0.05, y = 0.75),
                            "topright" = list(x = 0.75, y = 0.75),
                            "bottomleft" = list(x = 0.05, y = 0.05),
                            "bottomright" = list(x = 0.75, y = 0.05),
                            list(x = 0.75, y = 0.05))
    p <- ggdraw(p) +
      draw_grob(
        legend_grob,
        x = legend_coords$x,
        y = legend_coords$y,
        width = 0.2,
        height = 0.2
      )
  }
  return(p)
}
