#' Create a pathway glyph for one pathway
#'
#' Build a compact per-pathway summary graphic that can be embedded inside a
#' pathway graph or shown on its own. When more than one condition is supplied,
#' the glyph is typically used as the node image in [plot_subgraph()] or the
#' Shiny app. With `lollipop_plot = TRUE`, the glyph shows one lollipop per
#' condition. With `lollipop_plot = FALSE`, it shows a small bar chart instead.
#'
#' @details
#' `enrichment_scores` must be a data frame whose row names are pathway names and
#' whose columns are condition names. `pathway` must match one of those row
#' names. If `significance_scores` is supplied, it must use the same row/column
#' convention. The function returns both the in-memory magick image and the
#' temporary PNG path written to disk.
#'
#' This is a low-level helper that is useful when you want to inspect or reuse a
#' single glyph outside the full graph plotting workflow.
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_void scale_x_discrete geom_hline theme element_blank element_rect unit ylim geom_col scale_color_manual geom_point geom_segment scale_y_continuous scale_shape_manual scale_fill_manual expansion
#' @importFrom ragg agg_png
#' @importFrom magick image_read
#' @importFrom colorspace scale_fill_continuous_diverging scale_color_continuous_diverging
#' @importFrom grDevices dev.off
#' @param pathway Character scalar giving the pathway to draw. This must match a
#'   row name in `enrichment_scores`.
#' @param conditions Character vector of condition columns to include. If `NULL`
#'   or empty, all columns in `enrichment_scores` are used.
#' @param enrichment_scores Data frame of enrichment scores with pathways in rows
#'   and conditions in columns.
#' @param enrichment_limits Numeric vector of length 2 giving the y-axis/color
#'   scale limits. When `NULL`, limits are computed from the supplied scores.
#' @param significance_scores Optional data frame of significance values for the
#'   same pathways and conditions as `enrichment_scores`. Values less than or
#'   equal to `significance_cutoff` are marked as significant when
#'   `lollipop_significance = TRUE`.
#' @param glyph_size Numeric vector of length 2 giving glyph width and height in
#'   pixels.
#' @param res Numeric output resolution in dpi.
#' @param lollipop_plot Logical; if `TRUE`, use a lollipop glyph, otherwise use a
#'   small bar chart.
#' @param lollipop_colors Optional named character vector mapping condition names
#'   to colors. When `NULL`, colors are mapped to enrichment values.
#' @param lollipop_significance Logical; if `TRUE` and `significance_scores` is
#'   supplied, point shapes indicate significance.
#' @param significance_cutoff Numeric cutoff used to classify significance.
#' @param lollipop_size_scale Numeric scalar that scales both the lollipop stem
#'   width and point size. Values above `1` enlarge the glyph internals and
#'   values below `1` shrink them.
#' @return A list of length 2. The first element is a `magick-image` object and
#'   the second element is the temporary PNG filename used to create it.
#' @examples
#' pathway <- "Pathway1"
#' conditions <- c("Condition1", "Condition2")
#' enrichment_scores <- data.frame(
#'   Condition1 = c(1.2, -0.5),
#'   Condition2 = c(0.8, 1.1),
#'   row.names = c("Pathway1", "Pathway2")
#' )
#' significance_scores <- data.frame(
#'   Condition1 = c(0.01, 0.20),
#'   Condition2 = c(0.04, 0.30),
#'   row.names = c("Pathway1", "Pathway2")
#' )
#' glyph <- create_glyph_on_the_fly(
#'   pathway = pathway,
#'   conditions = conditions,
#'   enrichment_scores = enrichment_scores,
#'   enrichment_limits = c(-2, 2),
#'   significance_scores = significance_scores,
#'   lollipop_significance = TRUE
#' )
#' glyph[[1]]
#' @seealso [generate_glyph_images_cached()], [plot_subgraph()]
#' @export
create_glyph_on_the_fly <- function(pathway,
                                    conditions,
                                    enrichment_scores,
                                    enrichment_limits,
                                    significance_scores = NULL,
                                    glyph_size = c(80, 60),
                                    res = 96,
                                    lollipop_plot = TRUE,
                                    lollipop_colors = NULL,
                                    lollipop_significance = FALSE,
                                    significance_cutoff = 0.05,
                                    lollipop_size_scale = 1) {
  conditions <- resolve_conditions(enrichment_scores, conditions)
  enrichment_limits <- normalize_enrichment_limits(
    enrichment_scores = enrichment_scores[, conditions, drop = FALSE],
    enrichment_limits = enrichment_limits,
    zero_anchor = TRUE,
    symmetric = TRUE,
    padding = 0.06
  )
  
  if (!pathway %in% rownames(enrichment_scores)) {
    warning(paste("Pathway not found in enrichment_scores:", pathway))
    return(NULL)
  }
  
  scores_glyph <- enrichment_scores[pathway, conditions, drop = FALSE]
  if (all(is.na(scores_glyph))) {
    warning(paste("All enrichment scores are NA for pathway:", pathway))
    return(NULL)
  }
  
  df <- data.frame(
    Condition = factor(conditions, levels = conditions),
    Enrichment = as.numeric(scores_glyph[1, ]),
    stringsAsFactors = FALSE
  )
  
  show_significance <- isTRUE(lollipop_significance) && !is.null(significance_scores)
  if (show_significance) {
    significance_scores <- significance_scores[, conditions, drop = FALSE]
    sig_vals <- as.numeric(significance_scores[1, ])
    
    df$Significance <- ifelse(
      is.na(sig_vals),
      "NA",
      ifelse(
        sig_vals <= significance_cutoff & df$Enrichment >= 0,
        "S_up",
        ifelse(
          sig_vals <= significance_cutoff & df$Enrichment < 0,
          "S_down",
          "NS"
        )
      )
    )
    
    df$Significance <- factor(
      df$Significance,
      levels = c("S_up", "S_down", "NS", "NA")
    )
  }
  
  n_conditions <- max(length(conditions), 1L)
  res_ratio <- 96 / res
  
  lollipop_size_scale <- as.numeric(lollipop_size_scale)[1]
  if (!is.finite(lollipop_size_scale) || lollipop_size_scale <= 0) {
    lollipop_size_scale <- 1
  }
  
  seg_size <- max(
    (glyph_size[1] / 8 / n_conditions) * res_ratio * lollipop_size_scale,
    0.2
  )
  pt_size <- max(
    (glyph_size[1] / 6 / n_conditions) * res_ratio * lollipop_size_scale,
    0.5
  )
  point_stroke <- max(0.35 * sqrt(lollipop_size_scale), 0.25)
  hline_size <- max(0.5, 1.2 * res_ratio)
  use_custom_colors <- !is.null(lollipop_colors)
  
  if (lollipop_plot) {
    p_plot <- {
      if (use_custom_colors) {
        ggplot(df, aes(x = Condition, y = Enrichment, color = Condition))
      } else {
        ggplot(df, aes(x = Condition, y = Enrichment, color = Enrichment))
      }
    } +
      geom_hline(yintercept = 0, color = "grey45", linewidth = hline_size) +
      geom_segment(
        aes(xend = Condition, y = 0, yend = Enrichment),
        linewidth = seg_size,
        lineend = "round",
        na.rm = TRUE
      )
    
    if (show_significance) {
      if (use_custom_colors) {
        p_plot <- p_plot +
          geom_point(
            aes(shape = Significance, fill = Condition),
            color = "black",
            stroke = point_stroke,
            size = pt_size,
            na.rm = TRUE
          ) +
          scale_fill_manual(values = lollipop_colors, drop = FALSE)
      } else {
        p_plot <- p_plot +
          geom_point(
            aes(shape = Significance, fill = Enrichment),
            color = "black",
            stroke = point_stroke,
            size = pt_size,
            na.rm = TRUE
          ) +
          scale_fill_continuous_diverging(
            palette = "Berlin",
            limits = enrichment_limits
          )
      }
      
      p_plot <- p_plot +
        scale_shape_manual(
          values = c(
            "S_up" = 24,
            "S_down" = 25,
            "NS" = 21,
            "NA" = 4
          ),
          drop = FALSE
        )
    } else {
      p_plot <- p_plot + geom_point(size = pt_size, na.rm = TRUE)
    }
    
    p_plot <- p_plot +
      {
        if (use_custom_colors) {
          scale_color_manual(values = lollipop_colors, drop = FALSE)
        } else {
          scale_color_continuous_diverging(
            palette = "Berlin",
            limits = enrichment_limits
          )
        }
      } +
      ggplot2::coord_cartesian(ylim = enrichment_limits, clip = "off") +
      theme_void() +
      theme(
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0, 0, 0, 0), "pt")
      )
  } else {
    p_plot <- ggplot(df, aes(x = Condition, y = Enrichment, fill = Enrichment)) +
      geom_col(na.rm = TRUE, width = 0.92) +
      geom_hline(yintercept = 0, color = "grey45", linewidth = 0.5) +
      scale_x_discrete(expand = expansion(0, 0)) +
      scale_fill_continuous_diverging(palette = "Berlin", limits = enrichment_limits) +
      ggplot2::coord_cartesian(ylim = enrichment_limits, clip = "off") +
      theme_void() +
      theme(
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0, 0, 0, 0), "pt")
      )
  }
  
  img_file <- tempfile(fileext = ".png")
  agg_png(
    filename = img_file,
    width = glyph_size[1],
    height = glyph_size[2],
    units = "px",
    res = res,
    background = NA
  )
  print(p_plot)
  dev.off()
  
  img <- image_read(img_file)
  list(img, img_file)
}


#' Generate a legend for pathway glyphs
#'
#' Create a standalone legend that explains the glyph encoding used by
#' [create_glyph_on_the_fly()] and [plot_subgraph()]. This is primarily useful
#' when assembling publication figures or custom layouts outside the built-in
#' plotting helpers.
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal theme element_text element_blank element_rect unit geom_hline scale_y_continuous
#' @importFrom colorspace scale_fill_continuous_diverging scale_color_continuous_diverging
#' @param conditions Character vector of conditions shown in the glyphs.
#' @param enrichment_limits Numeric vector of length 2 defining the enrichment
#'   range represented by the legend.
#' @param reference Optional character scalar naming the reference condition used
#'   to center enrichment values. When supplied, the y-axis label states that the
#'   legend is relative to that reference.
#' @param lollipop_plot Logical; if `TRUE`, draw a lollipop-style legend,
#'   otherwise draw a bar-style legend.
#' @param lollipop_colors Optional named character vector mapping condition names
#'   to colors. When `NULL`, colors are mapped to enrichment values instead.
#' @param lollipop_significance Logical; if `TRUE`, show the significance shape
#'   encoding in the legend using a triangle for significant and a circle for
#'   non-significant values.
#' @param lollipop_size_scale Numeric scalar that scales both the lollipop stem
#'   width and point size in the legend.
#' @return A `ggplot` object that can be printed directly or combined with other
#'   plots.
#' @examples
#' generate_legend_plot(
#'   conditions = c("A", "B", "C"),
#'   enrichment_limits = c(-2, 2),
#'   lollipop_plot = TRUE
#' )
#' @seealso [create_glyph_on_the_fly()], [generate_relation_legend()]
#' @export
generate_legend_plot <- function(conditions,
                                 enrichment_limits,
                                 reference = NULL,
                                 lollipop_plot = TRUE,
                                 lollipop_colors = NULL,
                                 lollipop_significance = FALSE,
                                 lollipop_size_scale = 1) {
  if (is.null(conditions) || length(conditions) == 0) {
    return(NULL)
  }
  
  enrichment_limits <- normalize_enrichment_limits(
    enrichment_limits = enrichment_limits,
    zero_anchor = TRUE,
    symmetric = TRUE,
    padding = 0.06
  )
  
  legend_df <- data.frame(
    Condition = factor(conditions, levels = conditions),
    Enrichment = seq(
      enrichment_limits[1],
      enrichment_limits[2],
      length.out = length(conditions)
    )
  )
  
  y_label <- if (is.null(reference) || identical(reference, "None")) {
    "Enrichment score"
  } else {
    paste0("Enrichment score\nrelative to\n", reference)
  }
  
  legend_theme <- theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, size = 9, hjust = 1, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 9, face = "bold"),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.ticks = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(5, 5, 5, 5), "pt")
    )
  
  use_custom_colors <- !is.null(lollipop_colors)
  n_cond <- max(length(conditions), 1L)
  
  lollipop_size_scale <- as.numeric(lollipop_size_scale)[1]
  if (!is.finite(lollipop_size_scale) || lollipop_size_scale <= 0) {
    lollipop_size_scale <- 1
  }
  
  legend_seg_size <- max((9 / n_cond) * lollipop_size_scale, 0.3)
  legend_pt_size <- max((13 / n_cond) * lollipop_size_scale, 0.75)
  legend_stroke <- max(0.35 * sqrt(lollipop_size_scale), 0.25)
  
  if (lollipop_plot) {
    if (isTRUE(lollipop_significance)) {
      legend_df$Significance <- ifelse(
        legend_df$Enrichment > 0,
        "S_up",
        ifelse(legend_df$Enrichment < 0, "S_down", "NS")
      )
      
      legend_df$Significance <- factor(
        legend_df$Significance,
        levels = c("S_up", "S_down", "NS", "NA")
      )
    }
    
    legend_plot <- {
      if (use_custom_colors) {
        ggplot(legend_df, aes(x = Condition, y = Enrichment, color = Condition))
      } else {
        ggplot(legend_df, aes(x = Condition, y = Enrichment, color = Enrichment))
      }
    } +
      geom_hline(yintercept = 0, color = "grey45", linewidth = 0.8) +
      geom_segment(
        aes(xend = Condition, y = 0, yend = Enrichment),
        linewidth = legend_seg_size
      )
    
    if (isTRUE(lollipop_significance)) {
      if (use_custom_colors) {
        legend_plot <- legend_plot +
          geom_point(
            aes(shape = Significance, fill = Condition),
            size = legend_pt_size,
            stroke = legend_stroke,
            color = "black"
          ) +
          scale_fill_manual(values = lollipop_colors, drop = FALSE)
      } else {
        legend_plot <- legend_plot +
          geom_point(
            aes(shape = Significance, fill = Enrichment),
            size = legend_pt_size,
            stroke = legend_stroke,
            color = "black"
          ) +
          scale_fill_continuous_diverging(
            palette = "Berlin",
            limits = enrichment_limits
          )
      }
      
      legend_plot <- legend_plot +
        scale_shape_manual(
          values = c(
            "S_up" = 24,
            "S_down" = 25,
            "NS" = 21,
            "NA" = 4
          ),
          drop = FALSE
        )
    } else {
      legend_plot <- legend_plot + geom_point(size = legend_pt_size)
    }
    
    legend_plot <- legend_plot +
      scale_x_discrete(expand = expansion(0.1, 0.2)) +
      scale_y_continuous(expand = expansion(0.12, 0.2)) +
      {
        if (use_custom_colors) {
          scale_color_manual(
            values = lollipop_colors,
            name = "Condition",
            drop = FALSE
          )
        } else {
          scale_color_continuous_diverging(
            palette = "Berlin",
            limits = enrichment_limits,
            name = "Enrichment"
          )
        }
      } +
      ylab(y_label) +
      ggtitle("Glyph legend") +
      legend_theme
  } else {
    legend_plot <- ggplot(
      legend_df,
      aes(x = Condition, y = Enrichment, fill = Enrichment)
    ) +
      geom_bar(stat = "identity", width = 0.92) +
      geom_hline(yintercept = 0, color = "grey45", linewidth = 0.5) +
      scale_x_discrete(expand = expansion(0, 0)) +
      scale_fill_continuous_diverging(
        palette = "Berlin",
        limits = enrichment_limits,
        name = "Enrichment"
      ) +
      ylab(y_label) +
      ggtitle("Glyph legend") +
      legend_theme
  }
  
  legend_plot
}
#' Generate a legend for pathway relations
#'
#' Create a small legend describing the edge relation styles used in pathway
#' graphs. This is most helpful when you keep the relation-specific line colors
#' or linetypes in [plot_subgraph()] and want a matching standalone legend.
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_text theme_minimal theme element_blank element_text unit coord_cartesian
#' @param relations_present Character vector of relation labels to include, for
#'   example `c("is_a", "part_of")`.
#' @param relation_colors Named character vector mapping relation labels to edge
#'   colors.
#' @param relation_linetypes Named character vector mapping relation labels to
#'   `ggplot2` linetype values.
#' @param textsize Numeric text size for the relation labels.
#' @return A `ggplot` object containing the relation legend.
#' @examples
#' generate_relation_legend(c("is_a", "part_of", "regulates"))
#' @seealso [plot_subgraph()]
#' @export
generate_relation_legend <- function(relations_present,
                                     relation_colors = c(
                                       "is_a" = "black",
                                       "part_of" = "darkcyan",
                                       "regulates" = "orange",
                                       "positively_regulates" = "steelblue",
                                       "negatively_regulates" = "salmon"
                                     ),
                                     relation_linetypes = c(
                                       "is_a" = "solid",
                                       "part_of" = "dashed",
                                       "regulates" = "dotted",
                                       "positively_regulates" = "dotdash",
                                       "negatively_regulates" = "longdash"
                                     ),
                                     textsize = 5) {
  df <- data.frame(
    relation = relations_present,
    x = 1,
    y = rev(seq_along(relations_present)),
    color = relation_colors[relations_present],
    linetype = relation_linetypes[relations_present],
    stringsAsFactors = FALSE
  )
  
  ggplot(df, aes(x = x, y = y)) +
    geom_segment(
      aes(x = x - 0.4, xend = x + 0.4, yend = y,
          color = I(color), linetype = I(linetype)),
      linewidth = 1
    ) +
    geom_text(aes(x = x, y = y - 0.3, label = relation),
              color = "black", size = textsize, hjust = 0.5, fontface = "bold") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(5, 5, 5, 5), "pt")
    ) +
    coord_cartesian(clip = "off")
}

#' Generate cached glyph images for many pathways
#'
#' Create one glyph per pathway in a graph and store the rendered PNG files in a
#' cache directory. Existing glyphs are reused when the input data hash has not
#' changed, which makes repeated plotting much faster in iterative workflows and
#' in the Shiny app.
#'
#' @details
#' The return value is a named list of base64-encoded image URIs, suitable for
#' use with packages such as `ggimage` or `visNetwork`. Names correspond to
#' pathway labels in `V(g)$label`.
#'
#' `gene_sets` is currently accepted for interface consistency with the rest of
#' the package, but glyph generation itself depends on pathway labels, selected
#' conditions, and score matrices.
#'
#' @importFrom base64enc dataURI
#' @importFrom digest digest
#' @importFrom magick image_write
#' @param g `igraph` object representing the pathway graph.
#' @param gene_sets Named list of pathway gene sets.
#' @param enrichment_scores Data frame with pathway names in rows and condition
#'   names in columns.
#' @param conditions Character vector of conditions to render. If `NULL` or
#'   empty, all columns in `enrichment_scores` are used.
#' @param mapping Data frame containing at least `processed_name` and
#'   `exact_source` columns.
#' @param enrichment_limits Optional numeric vector of length 2 used for the
#'   glyph color and y-axis scale.
#' @param significance_scores Optional data frame of significance values aligned
#'   to `enrichment_scores`.
#' @param glyph_size Numeric vector of length 2 giving glyph width and height in
#'   pixels.
#' @param res Numeric resolution in dpi.
#' @param cache_dir Directory used to store rendered PNG files and the input hash.
#' @param progress Optional Shiny `Progress` object for reporting progress.
#' @param force_regenerate Logical; if `TRUE`, rerender all glyphs even when the
#'   cache hash matches.
#' @param lollipop_plot Logical; if `TRUE`, glyphs are rendered as lollipop
#'   plots, otherwise as bar plots.
#' @param lollipop_colors Optional named character vector mapping conditions to
#'   colors.
#' @param lollipop_significance Logical; if `TRUE`, encode significance in point
#'   shapes when significance data are available.
#' @param significance_cutoff Numeric threshold for significance.
#' @param lollipop_size_scale Numeric scalar that scales both the lollipop stem
#'   width and point size in generated glyphs.
#' @return Named list of base64 data URIs, one per pathway label in the graph.
#' @examples
#' g <- igraph::graph_from_data_frame(
#'   data.frame(from = "A", to = c("B", "C")),
#'   directed = TRUE
#' )
#' igraph::V(g)$label <- c("PathwayA", "PathwayB", "PathwayC")
#' gene_sets <- list(PathwayA = c("G1", "G2"), PathwayB = "G1", PathwayC = "G3")
#' mapping <- data.frame(
#'   processed_name = c("PathwayA", "PathwayB", "PathwayC"),
#'   exact_source = c("A", "B", "C")
#' )
#' enrichment_scores <- data.frame(
#'   cond1 = c(1.0, -0.5, 0.2),
#'   cond2 = c(0.8, 0.1, -0.4),
#'   row.names = c("PathwayA", "PathwayB", "PathwayC")
#' )
#' glyphs <- generate_glyph_images_cached(
#'   g = g,
#'   gene_sets = gene_sets,
#'   enrichment_scores = enrichment_scores,
#'   conditions = c("cond1", "cond2"),
#'   mapping = mapping,
#'   cache_dir = tempfile("glyph_cache_")
#' )
#' names(glyphs)
#' @seealso [create_glyph_on_the_fly()], [get_cache_directory()]
#' @export
generate_glyph_images_cached <- function(
    g, gene_sets, enrichment_scores, conditions, mapping,
    enrichment_limits = NULL, significance_scores = NULL,
    glyph_size = c(80, 60), res = 96,
    cache_dir = file.path(tempdir(), "glyph_cache"),
    progress = NULL, force_regenerate = FALSE,
    lollipop_plot = TRUE, lollipop_colors = NULL,
    lollipop_significance = FALSE, significance_cutoff = 0.05,
    lollipop_size_scale = 1
) {
  cache_dir <- get_cache_directory(cache_dir)
  conditions <- resolve_conditions(enrichment_scores, conditions)
  enrichment_limits <- normalize_enrichment_limits(
    enrichment_scores = enrichment_scores[, conditions, drop = FALSE],
    enrichment_limits = enrichment_limits,
    zero_anchor = TRUE,
    symmetric = TRUE,
    padding = 0.06
  )
  
  data_hash <- digest(
    list(
      enrichment_scores = enrichment_scores[, conditions, drop = FALSE],
      conditions = conditions,
      significance_scores = significance_scores,
      enrichment_limits = enrichment_limits,
      glyph_size = glyph_size,
      res = res,
      lollipop_plot = lollipop_plot,
      lollipop_colors = lollipop_colors,
      lollipop_significance = lollipop_significance,
      significance_cutoff = significance_cutoff,
      lollipop_size_scale = lollipop_size_scale
    ),
    algo = "spookyhash"
  )
  hash_filepath <- file.path(cache_dir, "data_hash.txt")
  
  regenerate <- TRUE
  if (file.exists(hash_filepath)) {
    stored_hash <- readLines(hash_filepath, warn = FALSE)
    if (length(stored_hash) > 0 && identical(stored_hash, data_hash) && !force_regenerate) {
      regenerate <- FALSE
    }
  }
  
  if (regenerate) {
    png_files <- list.files(cache_dir, pattern = "\\.png$", full.names = TRUE)
    if (length(png_files) > 0) {
      file.remove(png_files)
    }
    writeLines(data_hash, con = hash_filepath)
  }
  
  pathway_nodes <- unique(setdiff(V(g)$label, c("Root", NA_character_)))
  total_pathways <- length(pathway_nodes)
  sanitized_names <- vapply(pathway_nodes, sanitize_filename, FUN.VALUE = character(1))
  filepaths <- file.path(cache_dir, paste0(sanitized_names, ".png"))
  glyph_images <- vector("list", length = total_pathways)
  
  if (!is.null(progress) && total_pathways > 0) {
    progress$set(message = "Generating glyphs", value = 0)
    increment_value <- 1 / total_pathways
  }
  
  for (i in seq_len(total_pathways)) {
    pathway <- pathway_nodes[i]
    img_filepath <- filepaths[i]
    
    if (!file.exists(img_filepath) || regenerate) {
      pathway_significance <- NULL
      if (!is.null(significance_scores) && pathway %in% rownames(significance_scores)) {
        pathway_significance <- significance_scores[pathway, conditions, drop = FALSE]
      }
      
      img <- tryCatch({
        create_glyph_on_the_fly(
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
        )[[1]]
      }, error = function(e) {
        warning(paste("Error creating glyph for pathway:", pathway, "-", e$message))
        NULL
      })
      
      if (!is.null(img)) {
        image_write(img, path = img_filepath, format = "png")
        glyph_images[[i]] <- dataURI(file = img_filepath, mime = "image/png")
      } else {
        glyph_images[[i]] <- NA_character_
      }
    } else {
      glyph_images[[i]] <- dataURI(file = img_filepath, mime = "image/png")
    }
    
    if (!is.null(progress) && total_pathways > 0) {
      progress$inc(amount = increment_value, detail = paste("Processing:", pathway))
    }
  }
  
  if (!is.null(progress) && total_pathways > 0) {
    progress$close()
  }
  
  names(glyph_images) <- pathway_nodes
  glyph_images
}
