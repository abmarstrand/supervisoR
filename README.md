## supervisoR

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE.md) [![GitHub release](https://img.shields.io/github/v/release/abmarstrand/supervisoR?display_name=tag)](../../releases) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15590859.svg)](https://doi.org/10.5281/zenodo.15590859)

Visualize gene-set enrichment on **hierarchical pathway and ontology networks** in R. The package supports both **publication-oriented static figures** built with `ggplot2`/`ggraph` and **interactive exploration** with `Shiny` + `visNetwork`.

## What supervisoR is for

Most pathway-enrichment visualizations flatten results into ranked tables or dot plots. `supervisoR` keeps the **hierarchical structure** of the pathway database visible, so parent-child relationships remain interpretable while you compare many conditions at once.

The package is most useful when you already have:

1.  a table of pathway-level statistics, such as normalized enrichment scores or signed pathway scores,
2.  multiple contrasts, time points, or cohorts that you want to compare within the same pathway context.

## Main features

-   Hierarchy-aware subgraph plotting for pathway families or ontology branches.
-   Compact per-node glyphs for displaying many conditions simultaneously.
-   Single-condition plots with node color encoding.
-   Optional lollipop-style glyphs that improve separation between many conditions.
-   Optional edge-width encoding based on gene overlap between parent and child pathways.
-   Interactive browsing with search by pathway or gene in a Shiny app.
-   Built-in access to Reactome and GO biological process hierarchies through `msigdbr` plus included relation files.

## Installation

``` r
install.packages("remotes")
remotes::install_github("abmarstrand/supervisoR")

# Optional packages used in the real-data vignette
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "piano"))
```

## Quick start

``` r
library(supervisoR)

set.seed(42)
database <- "reactome"
species <- "Homo sapiens"

default_data <- load_default_data(database = database, species = species)

# Example multi-condition enrichment table
example_scores <- data.frame(
  Condition_A = rnorm(length(default_data$gene_sets), sd = 0.35),
  Condition_B = rnorm(length(default_data$gene_sets), sd = 0.35),
  Condition_C = rnorm(length(default_data$gene_sets), sd = 0.35)
)
rownames(example_scores) <- names(default_data$gene_sets)

p <- plot_subgraph(
  parent_geneset_name = "Metabolism Of Lipids",
  enrichment_scores = example_scores,
  conditions = colnames(example_scores),
  pathway_column = NULL,
  database = database,
  species = species,
  layout = "kk",
  hide_nodes_without_enrichment = TRUE,
  glyph_size = c(180, 130),
  lollipop_plot = TRUE,
  lollipop_colors = c(
    Condition_A = "#3B82F6",
    Condition_B = "#F59E0B",
    Condition_C = "#10B981"
  )
)

print(p)
```

## Inputs expectations

### Enrichment table

`plot_subgraph()` and `run_pathway_shiny_app()` expect a data frame whose rows are pathways and whose columns are conditions.

-   If pathway names are already stored as row names, use `pathway_column = NULL`.
-   If pathway names are in a regular column, keep `pathway_column` at the appropriate column index or name.
-   Values should usually be **signed** pathway statistics so direction is visible in the plots.

### Significance table

If you want significance-aware lollipop glyphs in static plots, provide a second data frame with the same pathway/condition layout. Typical inputs are adjusted p-values.

### Hierarchy inputs

By default, `load_default_data()` supplies Reactome or GO structures. You can also use your own hierarchy with:

-   a named list of gene sets,
-   a parent-child edge table with `parent` and `child` columns,
-   and optionally a translation layer linking readable names to hierarchy IDs.

## Static visualizations

`plot_subgraph()` is the main plotting function.

Useful options include:

-   `max_depth` to limit the visible hierarchy depth.
-   `reference` to show each condition relative to a selected baseline.
-   `adjust_edge_thickness = TRUE` to encode parent-child gene overlap.
-   `edge_percentage_labels = TRUE` to print overlap percentages on edges.
-   `repel_labels = TRUE` to improve label placement in dense figures.
-   `lollipop_plot = TRUE` for clearer multi-condition node glyphs.

## Interactive exploration

Use `run_pathway_shiny_app()` when you want to browse the full hierarchy instead of starting from a single parent node.

``` r
run_pathway_shiny_app(
  enrichment_scores = example_scores,
  conditions = colnames(example_scores),
  pathway_column = NULL,
  database = "reactome",
  species = "Homo sapiens",
  glyph_size = c(80, 60),
  res = 96
)
```

The app supports pathway search, gene search, drill-down into subgraphs, reset/back navigation, and automatic glyph caching.

## Custom hierarchies

If you already have your own ontology or pathway hierarchy, preprocess it once and then reuse the returned objects:

``` r
custom <- load_and_preprocess_gene_sets(
  gene_sets = my_gene_sets,
  pathways_relation = my_relations,
  translation_layer = my_translation
)

plot_subgraph(
  parent_geneset_name = "My Parent Pathway",
  enrichment_scores = my_scores,
  conditions = colnames(my_scores),
  pathway_column = NULL,
  mapping = custom$mapping,
  g = custom$pathway_graph,
  gene_sets = custom$gene_sets
)
```

## Documentation

The package currently ships with two vignettes.

-   `vignette("Introduction", package = "supervisoR")` introduces the plotting API with simulated data.
-   `vignette("Data_processing_example", package = "supervisoR")` shows a full RNA-seq to pathway-visualization workflow.

## Getting help

Please open an [issue](../../issues) with a minimal reproducible example when possible. Including the pathway name, selected conditions, and whether you used built-in or custom hierarchy data makes debugging much faster.

## Citation

If you use `supervisoR` in your work, please cite the Zenodo DOI shown above.

## License

MIT — see [LICENSE.md](LICENSE.md).
