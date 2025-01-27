# supervisoR
An R package for the visualization of enrichment scores in hierarchical pathway networks.

To install the package, the following packages are needed:
```
BiocManager::install(c("msigdbr", "igraph", "ggraph", "ggplot2", "ggrepel", "stringr",
    "colorspace", "ggimage", "magick", "cowplot", "magrittr", "dplyr", "ragg",
    "base64enc", "digest", "progressr", "readr", "rlang", "shiny", "visNetwork",
    "shinyWidgets")
)
```
Then run the following to install supervisoR:
```
remotes::install_git("abmarstrand/supervisoR")
```

# Running supervisoR
Two vignettes are included, showing how supervisoR is run. The **Introduction.Rmd** vignette gives a basic overview of the package functionality including both graph generation and how to run the shiny app. Here we utilize the included pathways and synthetic data to show the functionality.

The **Data_processing_example.Rmd** vignette gives a real world example utilizing real, bulk RNA sequencing data from GEO (GSE142025) from the 2019 Wang Group paper (Fan et al., [2019](https://doi.org/10.2337/db19-0204)). This data contains human data of Diabetic Kidney Disease split into early and late stage disease alongside healthy human samples.
We utilize the included pathway files to perform a pathway analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [piano](https://bioconductor.org/packages/release/bioc/html/piano.html).


