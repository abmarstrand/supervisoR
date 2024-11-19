library(testthat)
library(supervisoR)

test_that("create_glyph_on_the_fly creates a valid image", {
  pathway <- "Pathway1"
  conditions <- c("Condition1", "Condition2")
  enrichment_scores <- data.frame(
    Condition1 = c(1.2, -0.5),
    Condition2 = c(0.8, 1.1),
    row.names = c("Pathway1", "Pathway2")
  )
  enrichment_limits <- c(-1, 2)

  img_list <- create_glyph_on_the_fly(pathway, conditions, enrichment_scores, enrichment_limits)

  img <- img_list[[1]]  # Access the image object
  expect_true(!is.null(img))
  expect_true(inherits(img, "magick-image"))
})
