#' Create a glyph for a pathway on-the-fly
#'
#' Generates a glyph for a single pathway by creating a bar plot of enrichment scores across specified conditions.
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_void scale_x_discrete geom_hline theme element_blank
#' @importFrom ragg agg_png
#' @importFrom magick image_read
#' @importFrom colorspace scale_fill_continuous_diverging
#' @importFrom rlang .data
#' @importFrom grDevices dev.off
#' @param pathway Character string representing the pathway name.
#' @param conditions Character vector of conditions to include in the glyph.
#' @param enrichment_scores Data frame of enrichment scores (rows are pathways, columns are conditions).
#' @param enrichment_limits Numeric vector of length 2 specifying the enrichment limits.
#' @param glyph_size Numeric vector specifying the width and height of the glyph in pixels. Default is \code{c(80, 60)}.
#' @param res Numeric value specifying the resolution of the image in dpi. Default is 96.
#' @return A magick image representing the glyph for the pathway as well as the local path to the image.
#' @examples
#' # Example usage:
#' \dontrun{
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
create_glyph_on_the_fly <- function(pathway, conditions, enrichment_scores, enrichment_limits,
                                    glyph_size = c(80, 60), res = 96) {
  # Get enrichment scores for the pathway and specified conditions
  scores_glyph <- enrichment_scores[pathway, conditions, drop = FALSE]

  if (all(is.na(scores_glyph))) {
    warning(paste("All enrichment scores are NA for pathway:", pathway))
    return(NULL)
  }

  # Convert to data frame for plotting
  df <- data.frame(
    Condition = factor(conditions, levels = conditions),
    Enrichment = as.numeric(scores_glyph[1, ])
  )

  # Create the bar plot with transparent background
  p_plot <- ggplot(df, aes(x = .data$Condition,
                           y = .data$Enrichment,
                           fill = .data$Enrichment)) +
    geom_bar(stat = "identity", na.rm = TRUE) +
    scale_x_discrete(expand = expansion(0, 0)) +
    ylim(enrichment_limits) +
    scale_fill_continuous_diverging(palette = "Berlin", limits = enrichment_limits) +
    theme_void() +
    geom_hline(yintercept = 0, color = "darkred", linewidth = 0.5) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      plot.margin = unit(c(0, 0, 0, 0), "pt")
    )

  # Render the plot to a temporary file and read it as an image
  img_file <- tempfile(fileext = ".png")
  ragg::agg_png(
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
  img <- magick::image_read(img_file)
  img_list <- list(img, img_file)
  return(img_list)
}


#' Create a enrichment score legend.
#'
#' Generates a legend of enrichment scores across specified conditions.
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal theme element_text element_blank element_rect unit
#' @importFrom colorspace  scale_fill_continuous_diverging
#' @param conditions Character vector of conditions to include in the glyph.
#' @param enrichment_limits Numeric vector of length 2 specifying the enrichment limits.
#' @param reference Optional string describing which conditions should be used as a reference.
#' @return A ggplot object representing the legend.
#' @export
generate_legend_plot <- function(conditions, enrichment_limits, reference = NULL) {
  # Create a data frame for the legend
  legend_df <- data.frame(
    Condition = factor(conditions, levels = conditions),
    Enrichment = seq(enrichment_limits[1], enrichment_limits[2], length.out = length(conditions))
  )

  # Generate the legend plot
  legend_plot <- ggplot(legend_df, aes(x = .data$Condition,
                                       y = .data$Enrichment,
                                       fill = .data$Enrichment)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, color = "darkred", linewidth = 0.5) +
    scale_x_discrete(expand = expansion(0, 0)) +
    scale_fill_continuous_diverging(
      palette = "Berlin",
      limits = enrichment_limits,
      name = "Enrichment"
    ) +
    ylab(ifelse(reference == "None" | is.null(reference), "Enrichment Score", paste0("Enrichment Score\nRelative to\n", reference))) +
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
      plot.margin = unit(c(0, 0, 0, 0), "pt")
    )

  return(legend_plot)
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
#' @importFrom purrr map
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
#'
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
  force_regenerate = FALSE
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
  data_hash <- digest::digest(list(enrichment_scores = enrichment_scores, conditions = conditions), algo = "md5")
  hash_filepath <- file.path(cache_dir, "data_hash.txt")

  # Check if the hash file exists
  if (file.exists(hash_filepath)) {
    # Read the stored hash
    stored_hash <- readLines(hash_filepath, warn = FALSE)

    if (length(stored_hash) > 0 && stored_hash == data_hash && !force_regenerate) {
      message("Enrichment scores unchanged. Using cached glyphs.")
      # Proceed to load glyphs from cache without regenerating
    } else {
      message("Enrichment scores changed or force_regenerate is TRUE. Regenerating glyphs.")
      # Invalidate cache by removing existing glyphs
      png_files <- list.files(cache_dir, pattern = "\\.png$", full.names = TRUE)
      if (length(png_files) > 0) {
        file.remove(png_files)
      }
      # Update the hash file
      writeLines(data_hash, con = hash_filepath)
    }
  } else {
    # Hash file doesn't exist, create it
    message("No existing cache found. Generating glyphs.")
    writeLines(data_hash, con = hash_filepath)
  }

  # Retrieve all pathways present in the graph
  pathway_nodes <- setdiff(V(g)$label, "Root")
  total_pathways <- length(pathway_nodes)

  # Progress bar setup
  if (!is.null(progress)) {
    progress$set(message = "Generating Glyphs", value = 0)
  }

  # Iterate over each pathway
  for (pathway in pathway_nodes) {
    # Define the file path for the glyph image
    img_filename <- paste0(sanitize_filename(pathway), ".png")
    img_filepath <- file.path(cache_dir, img_filename)

    # Check if the image already exists in the cache
    if (file.exists(img_filepath)) {
      # Load the image and convert to base64 URI
      encoded_img <- base64enc::dataURI(file = img_filepath, mime = "image/png")
      glyph_images[[pathway]] <- encoded_img
    } else {
      # Create glyph image
      img <- tryCatch({
        create_glyph_on_the_fly(
          pathway = pathway,
          conditions = conditions,
          enrichment_scores = enrichment_scores,
          enrichment_limits = enrichment_limits,
          glyph_size = glyph_size,
          res = res
        )[[1]]
      }, error = function(e) {
        warning(paste("Error creating glyph for pathway:", pathway, "-", e$message))
        return(NULL)
      })

      if (!is.null(img)) {
        # Save the image to the cache directory
        magick::image_write(img, path = img_filepath, format = "png")

        # Convert to base64 URI
        encoded_img <- base64enc::dataURI(file = img_filepath, mime = "image/png")
        glyph_images[[pathway]] <- encoded_img
      } else {
        glyph_images[[pathway]] <- NA
      }
    }

    # Update progress if provided
    if (!is.null(progress)) {
      progress$inc(1 / total_pathways, detail = paste("Processing:", pathway))
    }
  }

  if (!is.null(progress)) {
    progress$close()
  }

  return(glyph_images)
}
