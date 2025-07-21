#' Create a glyph for a pathway on-the-fly
#'
#' Generates a glyph for a single pathway by creating a bar plot of enrichment scores across specified conditions.
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_void scale_x_discrete geom_hline theme element_blank ylim geom_col scale_color_manual geom_point geom_segment scale_y_discrete scale_shape_manual scale_fill_manual
#' @importFrom ragg agg_png
#' @importFrom magick image_read
#' @importFrom colorspace scale_fill_continuous_diverging scale_color_continuous_diverging
#' @importFrom rlang .data
#' @importFrom grDevices dev.off
#' @importFrom dplyr mutate_all case_when
#' @importFrom magrittr %>%
#' @param pathway Character string representing the pathway name.
#' @param conditions Character vector of conditions to include in the glyph.
#' @param enrichment_scores Data frame of enrichment scores (rows are pathways, columns are conditions).
#' @param enrichment_limits Numeric vector of length 2 specifying the enrichment limits.
#' @param glyph_size Numeric vector specifying the width and height of the glyph in pixels. Default is \code{c(80, 60)}.
#' @param res Numeric value specifying the resolution of the image in dpi. Default is 96.
#' @param lollipop_plot Bool specifying whether we should create lollipop plots (\code{TRUE}) or not (\code{FALSE}). Lollipop plots make large numbers of comparisons easier to compare than barplots, especially when combined with lollipop_colors, with a slight loss in enrichment strength fidelity.
#' @param lollipop_colors Optional named list, specifies the colors of conditions in the lollipop plot. Names must match condititions exactly
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
                                    lollipop_plot = T,
                                    lollipop_colors = NULL,
                                    lollipop_significance = F,
                                    significance_cutoff = 0.05) {
  # Get enrichment scores for the pathway and specified conditions
  scores_glyph <- enrichment_scores[pathway, conditions, drop = FALSE]

  if (all(is.na(scores_glyph))) {
    warning(paste("All enrichment scores are NA for pathway:", pathway))
    return(NULL)
  }

  # Convert to data frame for plotting
  df <- data.frame(
    Condition = factor(conditions, levels = conditions),
    Enrichment = as.numeric(scores_glyph[1, ]),
    stringsAsFactors = FALSE
  )
  if (!is.null(significance_scores)) {
    significance_scores <- significance_scores %>%
      mutate_all(~ case_when(. > significance_cutoff ~ "NS", T ~ "S"))
    df["Significance"] <- factor(significance_scores, levels = c("NS", "S"))
  }
  
  if(lollipop_plot == T) {
    p_plot <- {if(!is.null(lollipop_colors)) ggplot(df, aes(color = Condition)) else ggplot(df, aes(color = Enrichment))} +
      geom_hline(yintercept = 0, color = "darkred", linewidth = 2*(96/res)) +
      geom_segment(aes(x = Condition, xend = Condition, y = 0, yend = Enrichment), size = (glyph_size[1]/8/length(conditions))*(96/res)) +
      {if (!is.null(significance_scores)) {
        {if(!is.null(lollipop_colors)) {
            geom_point(aes(x = Condition, y = Enrichment, shape = Significance, fill = Condition),
                     color = "black",
                     size = (glyph_size[1]/6/length(conditions))*(96/res))
          } else {
            geom_point(aes(x = Condition, y = Enrichment, shape = Significance, fill = Enrichment),
                       color = "black",
                       size = (glyph_size[1]/6/length(conditions))*(96/res))
          }
        }} else {
          geom_point(aes(x = Condition, y = Enrichment), size = (glyph_size[1]/6/length(conditions))*(96/res))
        }
      } +
      ylim(enrichment_limits[1]*1.1, enrichment_limits[2]*1.1) +
      {if (!is.null(lollipop_colors)) {
        scale_color_manual(values = lollipop_colors)
      } else {
        scale_color_continuous_diverging(palette = "Berlin", limits = enrichment_limits) 
      }} +
      {if (!is.null(significance_scores)) {
        scale_shape_manual(values = c(23, 21), breaks = c("S", "NS"))
      }} +
      {if (!is.null(significance_scores) && !is.null(lollipop_colors)) {
        scale_fill_manual(values = lollipop_colors)
      }} +
      theme_void() +
      theme(
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0, 0, 0, 0), "pt")
      )
  } else {
    # Create the bar plot with transparent background
    p_plot <- ggplot(df, aes(x = Condition,
                             y = Enrichment,
                             fill = Enrichment)) +
      geom_col(na.rm = TRUE) +
      scale_x_discrete(expand = expansion(0, 0)) +
      ylim(enrichment_limits[1], enrichment_limits[2]*1.01) +
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


  # Render the plot to a temporary file and read it as an image
  img_file <- tempfile(fileext = ".png")
  agg_png(
    filename = img_file,
    width = glyph_size[1],
    height = glyph_size[2],
    units = "px",
    res = res,
    background = NA  # Use NA for transparent background
  )
  print(p_plot)
  dev.off()

  # Read the image back into R
  img <- image_read(img_file)
  img_list <- list(img, img_file)
  return(img_list)
}


#' Create an Enrichment Score Legend
#'
#' Generates a legend of enrichment scores across specified conditions.
#' This plot displays enrichment scores for each condition on a bar scale
#' and indicates the enrichment limits and any reference condition if provided.
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal theme element_text element_blank element_rect unit geom_hline
#' @importFrom colorspace scale_fill_continuous_diverging
#' @param conditions Character vector of conditions to include in the glyph.
#' @param enrichment_limits Numeric vector of length 2 specifying the enrichment limits.
#' @param reference Optional string describing which conditions should be used as a reference.
#' @param lollipop_plot Bool specifying whether we should create lollipop plots (\code{TRUE}) or not (\code{FALSE}). Lollipop plots make large numbers of comparisons easier to compare than barplots, especially when combined with lollipop_colors, with a slight loss in enrichment strength fidelity.
#' @param lollipop_colors Optional named list, specifies the colors of conditions in the lollipop plot. Names must match condititions exactly
#' @return A ggplot object representing the enrichment legend.
#' @export
generate_legend_plot <- function(conditions,
                                 enrichment_limits,
                                 reference = NULL,
                                 lollipop_plot = T,
                                 lollipop_colors = NULL) {
  legend_df <- data.frame(
    Condition = factor(conditions, levels = conditions),
    Enrichment = seq(enrichment_limits[1],
                     enrichment_limits[2],
                     length.out = length(conditions))
  )
  
  if (lollipop_plot) {
    legend_plot <- {if(!is.null(lollipop_colors)) ggplot(legend_df, aes(color = Condition)) else ggplot(legend_df, aes(color = Enrichment))} +
      geom_hline(yintercept = 0, color = "darkred", linewidth = 2) +
      geom_segment(aes(x = Condition, xend = Condition, y = 0, yend = Enrichment), size = 10/length(conditions)) +
      geom_point(aes(x = Condition, y = Enrichment), size = 15/length(conditions)) +
      scale_x_discrete(expand = expansion(0.1, 0.2)) +
      scale_y_discrete(expand = expansion(0.1, 0.2)) +
      {if (!is.null(lollipop_colors)) {
        scale_color_manual(values = lollipop_colors,
                           name = "Condition") 
      } else {
        scale_color_continuous_diverging(
          palette = "Berlin",
          limits = enrichment_limits,
          name = "Enrichment"
        ) 
      }} +
      ylab(ifelse(is.null(reference) || reference == "None",
                  "Enrichment Score",
                  paste0("Enrichment Score\nRelative to\n", reference))) +
      ggtitle("Legend") +
      theme_minimal(base_size = 9) +
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
  } else {
    legend_plot <- ggplot(legend_df, aes(x = Condition,
                                         y = Enrichment,
                                         fill = Enrichment)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0, color = "darkred", linewidth = 0.5) +
      scale_x_discrete(expand = expansion(0, 0)) +
      scale_fill_continuous_diverging(
        palette = "Berlin",
        limits = enrichment_limits,
        name = "Enrichment"
      ) +
      ylab(ifelse(is.null(reference) || reference == "None",
                  "Enrichment Score",
                  paste0("Enrichment Score\nRelative to\n", reference))) +
      ggtitle("Legend") +
      theme_minimal(base_size = 9) +
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
  }

  return(legend_plot)
}

#' Generate Relation Legend Plot
#'
#' Creates a standalone legend plot depicting the line types and colors associated
#' with different pathway relations. This plot can be combined with the enrichment
#' legend plot using \code{cowplot::plot_grid()} to display a unified legend.
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_text theme_minimal theme element_blank element_text unit coord_cartesian
#' @param relations_present A character vector of relations that are present in the current dataset (excluding NA).
#' @param relation_colors A named character vector mapping each relation to a color.
#' @param relation_linetypes A named character vector mapping each relation to a valid ggplot2 linetype
#'   (e.g., "solid", "dashed", "dotted", "dotdash", "longdash").
#' @param textsize Numeric value specifying the size of the relation labels. Default is 5.
#' @return A ggplot object representing the relation legend. Each relation is shown with its corresponding line style and color.
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
  # Create a data frame for plotting the relations
  df <- data.frame(
    relation = relations_present,
    x = 1,
    y = seq_along(relations_present),
    color = relation_colors[relations_present],
    linetype = relation_linetypes[relations_present],
    stringsAsFactors = FALSE
  )

  # Reverse the order so the first relation appears at the top
  df$y <- rev(df$y)

  rel_plot <- ggplot(df, aes(x = x, y = y)) +
    geom_segment(aes(x = x - 0.4, xend = x + 0.4, y = y, yend = y,
                     color = I(color), linetype = I(linetype)), size = 1) +
    geom_text(aes(x = x, y = y - 0.3, label = relation),
              color = "black", size = textsize, hjust = 0.5,
              fontface = "bold") +
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

  (rel_plot)
}


#' Generate Glyph Images with Caching
#'
#' Generates glyph images for each pathway and caches them to improve performance.
#' If the enrichment scores and conditions remain unchanged, cached images are used.
#'
#' @importFrom base64enc dataURI
#' @importFrom digest digest
#' @importFrom colorspace scale_fill_continuous_diverging
#' @importFrom magick image_write
#' @importFrom dplyr if_else
#' @importFrom progressr with_progress
#' @param g An \code{igraph} object representing the pathway graph.
#' @param gene_sets A named list where each element corresponds to a pathway and contains the associated genes.
#' @param enrichment_scores A \code{data.frame} with pathways as rows and conditions as columns, containing enrichment scores.
#' @param conditions A character vector specifying the conditions to visualize.
#' @param mapping A \code{data.frame} containing mapping information between pathways and their exact_source IDs.
#' @param enrichment_limits Optional numeric vector of length two specifying the minimum and maximum enrichment scores for visualization.
#' @param glyph_size Optional numeric vector specifying the width and height of glyphs in pixels. Default is \code{c(80, 60)}.
#' @param res Optional numeric value specifying the resolution of glyph images in DPI. Default is 96.
#' @param cache_dir Character string specifying the directory to store cached glyphs. Default is \code{"glyph_cache"}.
#' @param progress Optional Shiny progress object for updating progress bars.
#' @param force_regenerate Logical indicating whether to force regeneration of glyph images, ignoring the cache. Default is \code{FALSE}.
#' @param lollipop_plot Bool specifying whether we should create lollipop plots (\code{TRUE}) or not (\code{FALSE}). Lollipop plots make large numbers of comparisons easier to compare than barplots, especially when combined with lollipop_colors, with a slight loss in enrichment strength fidelity.
#' @param lollipop_colors Optional named list, specifies the colors of conditions in the lollipop plot. Names must match condititions exactly
#' @return A named list where each name corresponds to a pathway and contains the base64-encoded image URI.
#' @export
generate_glyph_images_cached <- function(
  g,
  gene_sets,
  enrichment_scores,
  conditions,
  mapping,
  enrichment_limits = NULL,
  glyph_size = c(80, 60),
  res = 96,
  cache_dir = file.path(tempdir(), "glyph_cache"),
  progress = NULL,
  force_regenerate = FALSE,
  lollipop_plot = T,
  lollipop_colors = NULL
) {
  # Initialize a list to store image URIs
  glyph_images <- list()

  # Ensure the cache directory exists
  cache_dir <- get_cache_directory(cache_dir)

  # Compute enrichment_limits if NULL
  if (is.null(enrichment_limits)) {
    enrichment_min <- min(enrichment_scores, na.rm = TRUE)
    enrichment_max <- max(enrichment_scores, na.rm = TRUE)
    enrichment_limits <- c(enrichment_min, enrichment_max)
    message(paste("Computed enrichment_limits:", enrichment_min, enrichment_max))
  }

  # Compute a hash of the enrichment_scores and conditions for cache invalidation
  data_hash <- digest(list(enrichment_scores = enrichment_scores, conditions = conditions), algo = "spookyhash")
  hash_filepath <- file.path(cache_dir, "data_hash.txt")

  # Check if we should regenerate glyphs
  regenerate <- TRUE
  # Check if the hash file exists
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
    # Invalidate old cache
    png_files <- list.files(cache_dir, pattern = "\\.png$", full.names = TRUE)
    if (length(png_files) > 0) {
      file.remove(png_files)
    }
    writeLines(data_hash, con = hash_filepath)
  }

  # Retrieve all pathways present in the graph
  pathway_nodes <- setdiff(V(g)$label, "Root")
  total_pathways <- length(pathway_nodes)
  
  sanitized_names <- vapply(pathway_nodes, sanitize_filename, FUN.VALUE = character(1))
  filepaths <- file.path(cache_dir, paste0(sanitized_names, ".png"))
  
  glyph_images <- vector("list", length = total_pathways)

  # Progress bar setup
  if (!is.null(progress)) {
    progress$set(message = "Generating Glyphs", value = 0)
    increment_value <- 1 / total_pathways
  }

  # Iterate over each pathway
  for (i in seq_len(total_pathways)) {
    # Define the file path for the glyph image
    pathway <- pathway_nodes[i]
    img_filepath <- filepaths[i]
    
    if (!file.exists(img_filepath) || regenerate) {
      img <- tryCatch({
        create_glyph_on_the_fly(
          pathway = pathway,
          conditions = conditions,
          enrichment_scores = enrichment_scores,
          enrichment_limits = enrichment_limits,
          glyph_size = glyph_size,
          res = res,
          lollipop_plot = lollipop_plot,
          lollipop_colors = lollipop_colors
        )[[1]]
      }, error = function(e) {
        warning(paste("Error creating glyph for pathway:", pathway, "-", e$message))
        return(NULL)
      })
      if (!is.null(img)) {
        image_write(img, path = img_filepath, format = "png")
        encoded_img <- base64enc::dataURI(file = img_filepath, mime = "image/png")
        glyph_images[[i]] <- encoded_img
      } else {
        glyph_images[[i]] <- NA
      }
    } else {
      # Use cached image
      encoded_img <- base64enc::dataURI(file = img_filepath, mime = "image/png")
      glyph_images[[i]] <- encoded_img
    }
    # Update progress if provided
    if (!is.null(progress)) {
      progress$inc(amount = increment_value, detail = paste("Processing:", pathway))
    }
  }
  if (!is.null(progress)) {
    progress$close()
  }
  names(glyph_images) <- pathway_nodes
  return(glyph_images)
}
