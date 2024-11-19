library(testthat)
library(supervisoR)
library(igraph)

test_that("get_child_pathways retrieves correct child pathways", {
  # Create a simple graph
  edges <- data.frame(
    from = c("A", "A", "B"),
    to = c("B", "C", "C"),
    stringsAsFactors = FALSE
  )
  g <- graph_from_data_frame(edges, directed = TRUE)

  # Define gene sets
  gene_sets <- list(
    "A" = c("Gene1", "Gene2"),
    "B" = c("Gene2", "Gene3"),
    "C" = c("Gene3", "Gene4")
  )

  # Define mapping
  mapping <- data.frame(
    processed_name = c("PathwayA", "PathwayB", "PathwayC"),
    exactSource = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )

  exactSource_to_name <- setNames(mapping$processed_name, mapping$exactSource)

  # Get child pathways of PathwayA
  children <- get_child_pathways("PathwayA", mapping, exactSource_to_name, g)

  expect_equal(sort(children), sort(c("PathwayB", "PathwayC")))

  # Get child pathways of PathwayB
  children_b <- get_child_pathways("PathwayB", mapping, exactSource_to_name, g)

  expect_equal(children_b, "PathwayC")

  # Get child pathways of PathwayC (should be none)
  children_c <- get_child_pathways("PathwayC", mapping, exactSource_to_name, g)

  expect_equal(children_c, character(0))
})
