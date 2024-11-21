#' Plot a Subgraph with Glyphs as Nodes
#'
#' Plots a subgraph with nodes represented by glyphs that visualize enrichment scores across specified conditions.
#' Glyphs are generated on-the-fly using the \code{create_glyph_on_the_fly} function, allowing dynamic comparisons
#' and efficient plotting even with a large number of conditions.
#'
#' @importFrom igraph subcomponent induced_subgraph V delete_vertices ends `V<-` E `E<-`
#' @importFrom ggraph ggraph create_layout geom_edge_link geom_node_label geom_node_text scale_edge_colour_continuous scale_edge_colour_viridis scale_edge_width geom_node_point
#' @importFrom ggplot2 theme_void scale_colour_manual ggplotGrob expansion ylim xlim ylab xlab ggtitle
#' @importFrom cowplot draw_grob ggdraw draw_grob
#' @importFrom ggimage geom_image
#' @importFrom magick image_write
#' @importFrom base64enc dataURI
#' @importFrom rlang abort
#' @importFrom grDevices png
#' @importFrom stats na.omit
#' @importFrom magrittr %>%
#' @importFrom colorspace scale_colour_continuous_diverging
#' @importFrom rlang .data
#' @param parent_geneset_name Name of the parent geneset.
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
#' @param color_by_depth Logical, whether to color edges by depth from the root node. Default is \code{FALSE}.
#' @param glyph_size Numeric vector specifying the width and height of the glyphs in pixels. Default is \code{c(80, 60)}.
#' @param res Numeric value specifying the resolution of the glyph images in dpi. Default is 96.
#' @param legend_position Optional character string ("topright", "topleft", "bottomleft" or "bottomright") determining final legend position.
#' @return A ggplot object representing the subgraph.
#' @export
plot_subgraph <- function(parent_geneset_name, enrichment_scores, conditions,
                          mapping, g, gene_sets = NULL,
                          layout = "kk", circular = FALSE, hide_nodes_without_enrichment = TRUE,
                          enrichment_limits = NULL, use_node_label = TRUE, reference = NULL,
                          adjust_edge_thickness = FALSE, edge_percentage_labels = FALSE, color_by_depth = FALSE,
                          glyph_size = c(80, 60), res = 96, legend_position = NULL) {
  # Validate inputs
  if (!parent_geneset_name %in% mapping$processed_name) {
    rlang::abort("Parent pathway not found in the mapping.")
  }

  parent_exact_source <- mapping$exact_source[match(parent_geneset_name, mapping$processed_name)]

  # Find the parent vertex in the graph
  parent_vertex <- which(V(g)$name == parent_exact_source)
  if (length(parent_vertex) == 0) {
    rlang::abort("Parent pathway not found in the graph.")
  }

  # Get all descendants (children) of the parent pathway
  descendants <- subcomponent(g, v = parent_vertex, mode = "out")

  # Induce subgraph
  sub_g <- induced_subgraph(g, vids = descendants)

  # Assign labels to subgraph vertices
  V(sub_g)$label <- mapping$processed_name[match(V(sub_g)$name, mapping$exact_source)]

  # Remove vertices with missing labels except the parent
  missing_labels <- is.na(V(sub_g)$label) & V(sub_g)$name != parent_exact_source
  if (any(missing_labels)) {
    sub_g <- delete_vertices(sub_g, V(sub_g)[missing_labels])
  }

  # Handle enrichment scores
  if (!is.null(enrichment_scores) && !is.null(conditions)) {
    # Adjust enrichment scores if reference is provided
    if (!is.null(reference) && reference %in% conditions) {
      enrichment_scores <- enrichment_scores - enrichment_scores[, reference, drop = TRUE]
      conditions <- setdiff(conditions, reference)
    }

    # Subset enrichment scores for the selected pathways and conditions
    scores <- enrichment_scores[V(sub_g)$label, conditions, drop = FALSE]

    # Handle missing scores
    missing_scores <- apply(scores, 1, function(x) all(is.na(x)))
    if (hide_nodes_without_enrichment) {
      sub_g <- induced_subgraph(sub_g, vids = V(sub_g)[!missing_scores])
      scores <- scores[!missing_scores, , drop = FALSE]
    }

    # Determine enrichment limits
    if (is.null(enrichment_limits)) {
      enrichment_limits <- range(scores, na.rm = TRUE)
    }
  }

  # If edge thickness or labels are to be adjusted, add overlap information
  if (adjust_edge_thickness || edge_percentage_labels) {
    sub_g <- add_overlap_to_edges(sub_g, gene_sets)
  }

  # If color_by_depth is TRUE, compute depth from the root for each edge
  if (color_by_depth) {
    distance_calc <- NULL
    E(sub_g)$distances <- 1
    for (edge_pw in ends(sub_g, E(sub_g))[, 2]) {
      distance_calc <- c(distance_calc, igraph::distances(g, parent_exact_source, edge_pw)[[1]])
    }
    E(sub_g)$distances <- distance_calc
    edge_colors <- scale_edge_colour_viridis(option = "G", aesthetics = "edge_colour",
                                             end = 0.75, guide = "none")
  } else {
    E(sub_g)$distances <- 1
    edge_colors <- scale_edge_colour_continuous(low = "gray50", high = "gray50",
                                                aesthetics = "edge_colour", guide = "none")
  }

  # Create the layout data with x and y
  layout_data <- create_layout(sub_g, layout = layout)

  # Determine if we have a single condition
  single_condition <- length(conditions) == 1

  if (!single_condition && !is.null(enrichment_scores) && !is.null(conditions)) {
    # Create glyph images on-the-fly
    glyphs <- list()
    for (pathway in V(sub_g)$label) {
      # Extract enrichment scores for the current pathway
      pathway_scores <- enrichment_scores[pathway, conditions, drop = FALSE]

      # Generate glyph image
      img <- create_glyph_on_the_fly(
        pathway = pathway,
        conditions = conditions,
        enrichment_scores = pathway_scores,
        enrichment_limits = enrichment_limits
      )[2]

      # Store the image
      if (!is.null(img)) {
        glyphs[[pathway]] <- img
      } else {
        glyphs[[pathway]] <- NA
      }
    }

    # Prepare node images
    node_images <- sapply(V(sub_g)$label, function(x) {
      glyphs[[x]][[1]]
    })

    # Add image column to layout_data
    layout_data$image <- node_images
  } else if (single_condition && !is.null(enrichment_scores)) {
    # For single condition, prepare enrichment scores for coloring
    single_condition_name <- conditions[1]
    layout_data$enrichment_score <- enrichment_scores[V(sub_g)$label, single_condition_name]
  }

  # Create the ggraph plot with layout_data
  p <- ggraph(layout_data) +
    geom_edge_link(
      aes(
        width = if (adjust_edge_thickness) .data$percent_overlap else 1,
        color = if (color_by_depth) .data$distances else 1
      ),
      show.legend = FALSE
    ) +
    scale_edge_width(range = c(0.5, 2)) +
    theme_void() +
    edge_colors

  if (!single_condition && !is.null(enrichment_scores) && !is.null(conditions)) {
    # Adjust glyph size based on number of nodes
    num_nodes <- nrow(layout_data)
    min_glyph_size <- 0.02  # Minimum glyph size
    max_glyph_size <- 0.1   # Maximum glyph size
    glyph_size_scale <- min(max_glyph_size, max(min_glyph_size, 0.5 / sqrt(num_nodes)))

    # Add glyphs as node images
    if (any(!is.na(layout_data$image))) {
      p <- p + ggimage::geom_image(aes(x = .data$x, y = .data$y, image = .data$image),
                                   size = glyph_size_scale, na.rm = TRUE)
    }

    if (any(is.na(layout_data$image))) {
      p <- p + geom_node_point(data = subset(layout_data, is.na(layout_data$image)),
                               aes(x = .data$x, y = .data$y),
                               color = "purple", size = glyph_size_scale * 25)
    }
  } else if (single_condition && !is.null(enrichment_scores)) {
    # Plot nodes as colored dots based on enrichment score
    p <- p + geom_node_point(aes(color = .data$enrichment_score), size = 5)

    # Define a continuous color scale for enrichment scores
    p <- p + scale_colour_continuous_diverging(palette = "Berlin", limits = enrichment_limits)
  }

  # Add node labels
  if (use_node_label) {
    p <- p + geom_node_label(aes(label = .data$label),
                             size = 3, repel = TRUE)
  } else {
    p <- p + geom_node_text(aes(label = .data$label),
                            size = 3, repel = TRUE)
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

  x_range <- x_max - x_min
  y_range <- y_max - y_min

  x_buffer <- x_range * 0.04  # Adjust buffer as needed
  y_buffer <- y_range * 0.04  # Adjust buffer as needed

  # Create the appropriate legend based on the condition
  if (!single_condition) {
    legend_plot <- generate_legend_plot(
      conditions = conditions,
      enrichment_limits = enrichment_limits,
      reference = reference
    )
    legend_grob <- ggplotGrob(legend_plot)

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
