library(testthat)
library(supervisoR)

test_that("load_default_data loads data correctly", {
  data <- load_default_data()

  expect_true(is.list(data))
  expect_true(all(c("gene_sets", "mapping", "exactSource_to_name", "pathway_graph", "relations") %in% names(data)))

  expect_true(is.list(data$gene_sets))
  expect_true(is.data.frame(data$mapping))
  expect_true(is.character(data$exactSource_to_name))
  expect_true(inherits(data$pathway_graph, "igraph"))
  expect_true(is.data.frame(data$relations))
})

test_that("load_and_preprocess_gene_sets works with and without translation_layer", {
  # Mock data
  gene_sets <- list(
    "Pathway1" = c("GeneA", "GeneB"),
    "Pathway2" = c("GeneC", "GeneD")
  )

  pathways_relation <- data.frame(
    parent = c("Parent1", "Parent1"),
    child = c("Pathway1", "Pathway2"),
    stringsAsFactors = FALSE
  )

  translation_layer <- data.frame(
    gene_sets_name = c("Pathway1", "Pathway2"),
    relation_name = c("Exact1", "Exact2"),
    stringsAsFactors = FALSE
  )

  # Without translation layer
  result_no_translation <- load_and_preprocess_gene_sets(gene_sets, pathways_relation)
  expect_true(is.list(result_no_translation))
  expect_true(all(c("gene_sets", "mapping", "exactSource_to_name", "pathway_graph") %in% names(result_no_translation)))

  # With translation layer
  result_with_translation <- load_and_preprocess_gene_sets(gene_sets, pathways_relation, translation_layer)
  expect_true(is.list(result_with_translation))
  expect_true(all(c("gene_sets", "mapping", "exactSource_to_name", "pathway_graph") %in% names(result_with_translation)))

  # Check mappings
  expect_equal(result_with_translation$mapping$exactSource, c("Exact1", "Exact2"))
})

