#' Create a glyph for a pathway on-the-fly
#'
#' Generates a glyph for a single pathway by creating a bar plot of enrichment scores across specified conditions.
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_void scale_x_discrete geom_hline theme element_blank element_rect unit ylim geom_col scale_color_manual geom_point geom_segment scale_y_continuous scale_shape_manual scale_fill_manual expansion
#' @importFrom ragg agg_png
#' @importFrom magick image_read
#' @importFrom colorspace scale_fill_continuous_diverging scale_color_continuous_diverging
#' @importFrom rlang .data
#' @importFrom grDevices dev.off
#' @importFrom dplyr mutate across everything case_when
#' @importFrom magrittr %>%
#' @param pathway Character string representing the pathway name.
#' @param conditions Character vector of conditions to include in the glyph.
#' @param enrichment_scores Data frame of enrichment scores (rows are pathways, columns are conditions).
#' @param enrichment_limits Numeric vector of length 2 specifying the enrichment limits.
#' @param significance_scores Optional data frame of significance scores (e.g. p-values). Must match the structure of enrichment_scores.
#' @param glyph_size Numeric vector specifying the width and height of the glyph in pixels. Default is \code{c(80, 60)}.
#' @param res Numeric value specifying the resolution of the image in dpi. Default is 96.
#' @param lollipop_plot Logical specifying whether we should create lollipop plots. Default is \code{TRUE}.
#' @param lollipop_colors Optional named list, specifies the colors of conditions in the lollipop plot.
#' @param lollipop_significance Logical specifying whether significance should be indicated. Default is \code{FALSE}.
#' @param significance_cutoff Numeric cutoff for significance. Default is 0.05.
#' @return A magick image representing the glyph for the pathway as well as the local path to the image.
#' @examples
#' # Example usage:
#' \donttest{
#' pathway <- "Pathway1"
#' conditions <- c("Condition1", "Condition2")
#' enrichment_scores <- data.frame(
#'   Condition1 = c(1.2, -0.5),
#'   Condition2 = c(0.8, 1.1)
#' )
#' rownames(enrichment_scores) <- c("Pathway1", "Pathway2")
#' enrichment_limits <- c(-1, 2)
#' glyph <- create_glyph_on_the_fly(pathway, conditions,
#'   enrichment_scores, enrichment_limits)
#' }
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
                                    significance_cutoff = 0.05) {
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
  
  n_conditions <- length(conditions)
  res_ratio <- 96 / res
  seg_size <- (glyph_size[1] / 8 / n_conditions) * res_ratio
  pt_size <- (glyph_size[1] / 6 / n_conditions) * res_ratio
  hline_size <- 2 * res_ratio
  ylim_lo <- enrichment_limits[1] * 1.1
  ylim_hi <- enrichment_limits[2] * 1.1
  use_custom_colors <- !is.null(lollipop_colors)
  
  if (!is.null(significance_scores)) {
    sig_vals <- as.data.frame(lapply(significance_scores, function(x) {
      ifelse(x > significance_cutoff, "NS", "S")
    }))
    df[["Significance"]] <- factor(unlist(sig_vals[1, ]), levels = c("NS", "S"))
  }
  
  if (lollipop_plot) {
    p_plot <- {
      if (use_custom_colors) {
        ggplot(df, aes(color = Condition))
      } else {
        ggplot(df, aes(color = Enrichment))
      }
    } +
      geom_hline(yintercept = 0, color = "darkred", linewidth = hline_size) +
      geom_segment(aes(x = Condition, xend = Condition, y = 0, yend = Enrichment),
                            linewidth = seg_size) +
      {
        if (!is.null(significance_scores)) {
          if (use_custom_colors) {
            geom_point(aes(x = Condition, y = Enrichment,
                                             shape = Significance, fill = Condition),
                                color = "black", size = pt_size)
          } else {
            geom_point(aes(x = Condition, y = Enrichment,
                                             shape = Significance, fill = Enrichment),
                                color = "black", size = pt_size)
          }
        } else {
          geom_point(aes(x = Condition, y = Enrichment), size = pt_size)
        }
      } +
      ylim(ylim_lo, ylim_hi) +
      {
        if (use_custom_colors) {
          scale_color_manual(values = lollipop_colors)
        } else {
          scale_color_continuous_diverging(palette = "Berlin", limits = enrichment_limits)
        }
      } +
      { if (!is.null(significance_scores)) scale_shape_manual(values = c(23, 21), breaks = c("S", "NS")) } +
      { if (!is.null(significance_scores) && use_custom_colors) scale_fill_manual(values = lollipop_colors) } +
      theme_void() +
      theme(
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0, 0, 0, 0), "pt")
      )
  } else {
    p_plot <- ggplot(df, aes(x = Condition, y = Enrichment, fill = Enrichment)) +
      geom_col(na.rm = TRUE) +
      scale_x_discrete(expand = expansion(0, 0)) +
      ylim(enrichment_limits[1], enrichment_limits[2] * 1.01) +
      scale_fill_continuous_diverging(palette = "Berlin", limits = enrichment_limits) +
      theme_void() +
      geom_hline(yintercept = 0, color = "darkred", linewidth = 0.5) +
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
  img_list <- list(img, img_file)
  return(img_list)
}


#' Create an Enrichment Score Legend
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal theme element_text element_blank element_rect unit geom_hline scale_y_continuous
#' @importFrom colorspace scale_fill_continuous_diverging
#' @param conditions Character vector of conditions to include in the glyph.
#' @param enrichment_limits Numeric vector of length 2 specifying the enrichment limits.
#' @param reference Optional string describing which conditions should be used as a reference.
#' @param lollipop_plot Logical specifying whether we should create lollipop plots.
#' @param lollipop_colors Optional named list, specifies the colors of conditions in the lollipop plot.
#' @return A ggplot object representing the enrichment legend.
#' @export
generate_legend_plot <- function(conditions,
                                 enrichment_limits,
                                 reference = NULL,
                                 lollipop_plot = TRUE,
                                 lollipop_colors = NULL) {
  legend_df <- data.frame(
    Condition = factor(conditions, levels = conditions),
    Enrichment = seq(enrichment_limits[1], enrichment_limits[2], length.out = length(conditions))
  )
  
  y_label <- ifelse(is.null(reference) || reference == "None",
                    "Enrichment Score",
                    paste0("Enrichment Score\nRelative to\n", reference))
  
  legend_theme <- theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, size = 10, hjust = 1, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(5, 5, 5, 5), "pt")
    )
  
  use_custom_colors <- !is.null(lollipop_colors)
  n_cond <- length(conditions)
  
  if (lollipop_plot) {
    legend_plot <- {
      if (use_custom_colors) ggplot(legend_df, aes(color = Condition))
      else ggplot(legend_df, aes(color = Enrichment))
    } +
      geom_hline(yintercept = 0, color = "darkred", linewidth = 2) +
      geom_segment(aes(x = Condition, xend = Condition, y = 0, yend = Enrichment),
                            linewidth = 10 / n_cond) +
      geom_point(aes(x = Condition, y = Enrichment), size = 15 / n_cond) +
      scale_x_discrete(expand = expansion(0.1, 0.2)) +
      scale_y_continuous(expand = expansion(0.1, 0.2)) +
      {
        if (use_custom_colors) scale_color_manual(values = lollipop_colors, name = "Condition")
        else scale_color_continuous_diverging(palette = "Berlin", limits = enrichment_limits, name = "Enrichment")
      } +
      ylab(y_label) +
      ggtitle("Legend") +
      legend_theme
  } else {
    legend_plot <- ggplot(legend_df, aes(x = Condition, y = Enrichment, fill = Enrichment)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0, color = "darkred", linewidth = 0.5) +
      scale_x_discrete(expand = expansion(0, 0)) +
      scale_fill_continuous_diverging(palette = "Berlin", limits = enrichment_limits, name = "Enrichment") +
      ylab(y_label) +
      ggtitle("Legend") +
      legend_theme
  }
  return(legend_plot)
}

#' Generate Relation Legend Plot
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_text theme_minimal theme element_blank element_text unit coord_cartesian
#' @param relations_present A character vector of relations that are present in the current dataset.
#' @param relation_colors A named character vector mapping each relation to a color.
#' @param relation_linetypes A named character vector mapping each relation to a ggplot2 linetype.
#' @param textsize Numeric value specifying the size of the relation labels. Default is 5.
#' @return A ggplot object representing the relation legend.
#' @export
generate_relation_legend <- function(relations_present,
                                     relation_colors = c("is_a" = "black", "part_of" = "darkcyan", "regulates" = "orange",
                                                         "positively_regulates" = "steelblue", "negatively_regulates" = "salmon"),
                                     relation_linetypes = c("is_a" = "solid", "part_of" = "dashed", "regulates" = "dotted",
                                                            "positively_regulates" = "dotdash", "negatively_regulates" = "longdash"),
                                     textsize = 5) {
  df <- data.frame(
    relation = relations_present, x = 1, y = seq_along(relations_present),
    color = relation_colors[relations_present], linetype = relation_linetypes[relations_present],
    stringsAsFactors = FALSE
  )
  df$y <- rev(df$y)
  
  ggplot(df, aes(x = x, y = y)) +
    geom_segment(aes(x = x - 0.4, xend = x + 0.4, y = y, yend = y,
                                       color = I(color), linetype = I(linetype)), linewidth = 1) +
    geom_text(aes(x = x, y = y - 0.3, label = relation),
                       color = "black", size = textsize, hjust = 0.5, fontface = "bold") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text = element_blank(), axis.ticks = element_blank(),
                   legend.position = "none", plot.margin = unit(c(5, 5, 5, 5), "pt")) +
    coord_cartesian(clip = "off")
}


#' Generate Glyph Images with Caching
#'
#' @importFrom base64enc dataURI
#' @importFrom digest digest
#' @importFrom magick image_write
#' @param g An igraph object representing the pathway graph.
#' @param gene_sets A named list of gene sets.
#' @param enrichment_scores A data.frame with pathways as rows and conditions as columns.
#' @param conditions A character vector specifying the conditions to visualize.
#' @param mapping A data.frame containing mapping information.
#' @param enrichment_limits Optional numeric vector of length two.
#' @param significance_scores Optional data frame of significance scores.
#' @param glyph_size Optional numeric vector for glyph dimensions. Default is c(80, 60).
#' @param res Optional numeric resolution in DPI. Default is 96.
#' @param cache_dir Character string for cache directory.
#' @param progress Optional Shiny progress object.
#' @param force_regenerate Logical to force regeneration. Default is FALSE.
#' @param lollipop_plot Logical for lollipop plot style. Default is TRUE.
#' @param lollipop_colors Optional named list for condition colors.
#' @param lollipop_significance Logical for significance display. Default is FALSE.
#' @param significance_cutoff Numeric cutoff. Default is 0.05.
#' @return A named list of base64-encoded image URIs.
#' @export
generate_glyph_images_cached <- function(
    g, gene_sets, enrichment_scores, conditions, mapping,
    enrichment_limits = NULL, significance_scores = NULL,
    glyph_size = c(80, 60), res = 96,
    cache_dir = file.path(tempdir(), "glyph_cache"),
    progress = NULL, force_regenerate = FALSE,
    lollipop_plot = TRUE, lollipop_colors = NULL,
    lollipop_significance = FALSE, significance_cutoff = 0.05
) {
  cache_dir <- get_cache_directory(cache_dir)
  
  if (is.null(enrichment_limits)) {
    enrichment_limits <- c(min(enrichment_scores, na.rm = TRUE), max(enrichment_scores, na.rm = TRUE))
    message(paste("Computed enrichment_limits:", enrichment_limits[1], enrichment_limits[2]))
  }
  
  data_hash <- digest(list(enrichment_scores = enrichment_scores,
                                   conditions = conditions,
                                   significance_scores = significance_scores),
                              algo = "spookyhash")
  hash_filepath <- file.path(cache_dir, "data_hash.txt")
  
  regenerate <- TRUE
  if (file.exists(hash_filepath)) {
    stored_hash <- readLines(hash_filepath, warn = FALSE)
    if (length(stored_hash) > 0 && stored_hash == data_hash && !force_regenerate) {
      message("Enrichment scores unchanged. Using cached glyphs.")
      regenerate <- FALSE
    } else {
      message("Enrichment scores changed or force_regenerate is TRUE. Regenerating glyphs.")
    }
  } else {
    message("No existing cache found. Generating glyphs.")
  }
  
  if (regenerate) {
    png_files <- list.files(cache_dir, pattern = "\\.png$", full.names = TRUE)
    if (length(png_files) > 0) file.remove(png_files)
    writeLines(data_hash, con = hash_filepath)
  }
  
  pathway_nodes <- setdiff(V(g)$label, "Root")
  total_pathways <- length(pathway_nodes)
  sanitized_names <- vapply(pathway_nodes, sanitize_filename, FUN.VALUE = character(1))
  filepaths <- file.path(cache_dir, paste0(sanitized_names, ".png"))
  glyph_images <- vector("list", length = total_pathways)
  
  if (!is.null(progress)) {
    progress$set(message = "Generating Glyphs", value = 0)
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
          pathway = pathway, conditions = conditions,
          enrichment_scores = enrichment_scores, enrichment_limits = enrichment_limits,
          significance_scores = pathway_significance, glyph_size = glyph_size, res = res,
          lollipop_plot = lollipop_plot, lollipop_colors = lollipop_colors,
          lollipop_significance = lollipop_significance, significance_cutoff = significance_cutoff
        )[[1]]
      }, error = function(e) {
        warning(paste("Error creating glyph for pathway:", pathway, "-", e$message))
        return(NULL)
      })
      if (!is.null(img)) {
        image_write(img, path = img_filepath, format = "png")
        glyph_images[[i]] <- dataURI(file = img_filepath, mime = "image/png")
      } else {
        glyph_images[[i]] <- NA
      }
    } else {
      glyph_images[[i]] <- dataURI(file = img_filepath, mime = "image/png")
    }
    if (!is.null(progress)) progress$inc(amount = increment_value, detail = paste("Processing:", pathway))
  }
  if (!is.null(progress)) progress$close()
  names(glyph_images) <- pathway_nodes
  return(glyph_images)
}