#' Sanitize filename
#'
#' Replaces invalid characters in a filename with underscores.
#'
#' @param name Character string representing the filename to sanitize.
#' @return A sanitized filename string.
#' @examples
#' sanitize_filename("path/to/file name with spaces.txt")
#' @export
sanitize_filename <- function(name) {
  name <- gsub("[^A-Za-z0-9_]", "_", name)
  return(name)
}

#' Check if a string is a file path
#'
#' Determines if a given string represents a file path by checking for path separators.
#'
#' @param name Character string to check.
#' @return Logical value indicating whether the string is a path.
#' @examples
#' is_path("path/to/file.txt") # Returns TRUE
#' is_path("filename") # Returns FALSE
#' @export
is_path <- function(name) {
  return(grepl("[/\\\\]", name))
}

#' Get or Create the Cache Directory
#'
#' Ensures that the specified cache directory exists. If it doesn't, the function creates it.
#'
#' @param cache_dir Character string specifying the directory to store cached glyphs. Default is \code{"glyph_cache"}.
#'
#' @return A character string representing the path to the cache directory.
#'
#' @examples
#' cache_path <- get_cache_directory("my_cache")
#'
#' @export
get_cache_directory <- function(cache_dir =  file.path(tempdir(), "glyph_cache")) {
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  }
  return(cache_dir)
}
