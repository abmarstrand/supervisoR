library(testthat)
library(supervisoR)
library(igraph)

test_that("add_overlap_to_edges adds correct overlap information", {
  # Create a simple graph
  edges <- data.frame(
    from = c("A", "A", "B"),
    to = c("B", "C", "C"),
    stringsAsFactors = FALSE
  )
  g <- igraph::graph_from_data_frame(edges, directed = TRUE)
  # Add label attribute to vertices
  igraph::V(g)$label <- c("PathwayA", "PathwayB", "PathwayC")
  # Define gene sets
  gene_sets <- list(
    "PathwayA" = c("Gene1", "Gene2"),
    "PathwayB" = c("Gene2", "Gene3"),
    "PathwayC" = c("Gene3", "Gene4")
  )
  V(g)$genes <- gene_sets[V(g)$label]
  # Add overlap information
  g <- add_overlap_to_edges(g, gene_sets)
  # Check overlap
  expected_overlaps <- c(1, 0, 1)  # Correct expectations
  expect_equal(igraph::E(g)$overlap, expected_overlaps)
  # Check percent_overlap
  expected_percent_overlap <- c(50, 0, 50)  # Correct expectations
  expect_equal(igraph::E(g)$percent_overlap, expected_percent_overlap)
})
