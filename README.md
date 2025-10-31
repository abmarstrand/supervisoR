## supervisoR
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE.md)
[![GitHub release](https://img.shields.io/github/v/release/abmarstrand/supervisoR?display_name=tag)](../../releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15590859.svg)](https://doi.org/10.5281/zenodo.15590859)

---
Visualize gene-set enrichment on **hierarchical pathway and ontology networks** in R — publishable static figures (ggplot2/ggraph) and interactive exploration via Shiny + visNetwork.
## Features
- **Hierarchy-aware plots:** Preserve biological structure by mapping enrichment onto pathway trees/DAGs.
- **Two modes:**  
  - *Static* figures for papers (ggplot2/ggraph)  
  - *Interactive* app for exploration (Shiny/visNetwork)
- **Scales to many conditions:** Layouts and encodings (incl. lollipop-style plots) to view many contrasts simultaneously.
---
## Installation
```r
# install supervisoR
install.packages("remotes")   # if needed
remotes::install_github("abmarstrand/supervisoR")

# optional: packages used in the real-data vignette
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "piano"))
```
---
## Quick start

```r
library(supervisoR)

# Explore the introductory vignette:
vignette("Introduction", package = "supervisoR")

# Or open the real-data workflow:
vignette("Data_processing_example", package = "supervisoR")
```

These vignettes walk through graph generation and how to run the Shiny app on both toy and real datasets.

---
## Inputs & Outputs

**Inputs (conceptual):**
1. A table of enrichment results per pathway × condition (e.g., Normalized Enrichment Scores from GSEA/fGSEA).
2. A pathway or ontology hierarchy (node and edge tables) to plot on.

**Outputs:**
- Publication-ready ggplot/ggraph figures.
- An interactive Shiny/visNetwork app for exploration and sharing.

---

## Documentation

- **Help pages:** `?supervisoR` and individual function help once loaded.
- **Vignettes:** `Introduction` and `Data_processing_example` cover end-to-end usage.

---

## Getting help

- Found a bug or have a feature request? Please open an [Issue](../../issues) with a minimal reproducible example if possible.

---

## Changelog

See [Releases](../../releases) for version history and highlights.

---

## Citing

If you use *supervisoR* in your work, please cite the Zenodo DOI:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15590859.svg)](https://doi.org/10.5281/zenodo.15590859)

---

## License

MIT — see [LICENSE.md](LICENSE.md).


