library(testthat)
library(supervisoR)
library(shinytest2)

test_that("Shiny app launches without errors", {
  # Mock data
  gene_sets <- list(
    "PathwayA" = c("Gene1", "Gene2"),
    "PathwayB" = c("Gene2", "Gene3"),
    "PathwayC" = c("Gene3", "Gene4")
  )

  pathways_relation <- data.frame(
    parent = c("A", "A"),
    child = c("B", "C"),
    stringsAsFactors = FALSE
  )

  mapping <- data.frame(
    processed_name = c("PathwayA", "PathwayB", "PathwayC"),
    exactSource = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )

  exactSource_to_name <- setNames(mapping$processed_name, mapping$exactSource)

  g <- igraph::graph_from_data_frame(pathways_relation, directed = TRUE)

  enrichment_scores <- data.frame(
    Condition1 = c(1.2, -0.5, 0.8),
    Condition2 = c(0.8, 1.1, -0.3),
    row.names = c("PathwayA", "PathwayB", "PathwayC")
  )

  # Launch the Shiny app in a separate process using shinytest2
  app <- shinytest2::AppDriver$new(
    app = run_pathway_shiny_app(
      enrichment_scores = enrichment_scores,
      conditions = c("Condition1", "Condition2"),
      mapping = mapping,
      exactSource_to_name = exactSource_to_name,
      g = g,
      gene_sets = gene_sets
    ),
    name = "supervisoR_shiny_app"
  )

  # Perform a simple interaction test
  expect_error(app$stop(), NA)  # Expect no error on stop
})
