#' Sanitize a string for use as a filename
#'
#' Replace characters outside letters, digits, underscores, and hyphens with
#' underscores. This is mainly used to create stable cache filenames for pathway
#' glyph images.
#'
#' @param name Character vector to sanitize.
#' @return Character vector with disallowed filename characters replaced by
#'   underscores.
#' @examples
#' sanitize_filename("path/to/file name with spaces.txt")
#' @seealso [get_cache_directory()]
#' @export
sanitize_filename <- function(name) {
  name <- gsub("[^A-Za-z0-9_\\-]", "_", name)
  return(name)
}

#' Test whether a string looks like a file path
#'
#' Check whether a string contains Unix or Windows path separators.
#'
#' @param name Character vector to test.
#' @return Logical vector; `TRUE` indicates that at least one path separator was
#'   detected.
#' @examples
#' is_path("path/to/file.txt")
#' is_path("filename")
#' @export
is_path <- function(name) {
  return(grepl("[/\\\\]", name))
}

#' Create or return a glyph cache directory
#'
#' Ensure that a directory exists and return its path. This is a small utility
#' used by [generate_glyph_images_cached()] and the Shiny app.
#'
#' @param cache_dir Directory path to use for glyph caching.
#' @return Character scalar giving the cache directory path.
#' @examples
#' cache_path <- get_cache_directory(file.path(tempdir(), "my_cache"))
#' dir.exists(cache_path)
#' @seealso [generate_glyph_images_cached()]
#' @export
get_cache_directory <- function(cache_dir = file.path(tempdir(), "glyph_cache")) {
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  }
  return(cache_dir)
}

# Internal helper: compute robust enrichment limits with optional zero anchoring.
normalize_enrichment_limits <- function(enrichment_scores = NULL,
                                        enrichment_limits = NULL,
                                        zero_anchor = TRUE,
                                        symmetric = FALSE,
                                        padding = 0.08) {
  if (is.null(enrichment_limits)) {
    if (is.null(enrichment_scores)) {
      stop("Either enrichment_scores or enrichment_limits must be supplied.")
    }
    values <- suppressWarnings(as.numeric(unlist(enrichment_scores, use.names = FALSE)))
    values <- values[is.finite(values)]
    if (length(values) == 0) {
      enrichment_limits <- c(-1, 1)
    } else {
      enrichment_limits <- range(values)
    }
  }
  
  enrichment_limits <- as.numeric(enrichment_limits)
  if (length(enrichment_limits) != 2 || any(!is.finite(enrichment_limits))) {
    enrichment_limits <- c(-1, 1)
  }
  
  if (diff(enrichment_limits) == 0) {
    delta <- max(abs(enrichment_limits[1]) * 0.05, 1e-6)
    enrichment_limits <- enrichment_limits + c(-delta, delta)
  }
  
  if (zero_anchor) {
    enrichment_limits[1] <- min(enrichment_limits[1], 0)
    enrichment_limits[2] <- max(enrichment_limits[2], 0)
  }
  
  if (symmetric) {
    max_abs <- max(abs(enrichment_limits))
    enrichment_limits <- c(-max_abs, max_abs)
  }
  
  span <- diff(enrichment_limits)
  pad <- if (isTRUE(padding) || is.numeric(padding)) span * padding else 0
  if (!is.finite(pad) || pad <= 0) {
    pad <- max(abs(enrichment_limits)) * 0.05
  }
  
  enrichment_limits + c(-pad, pad)
}

# Internal helper: return the first matching column names and fail early on missing inputs.
resolve_conditions <- function(enrichment_scores, conditions) {
  if (is.null(conditions) || length(conditions) == 0) {
    return(colnames(enrichment_scores))
  }
  
  missing_conditions <- setdiff(conditions, colnames(enrichment_scores))
  if (length(missing_conditions) > 0) {
    stop(
      "The following conditions were not found in enrichment_scores: ",
      paste(missing_conditions, collapse = ", ")
    )
  }
  
  conditions
}
