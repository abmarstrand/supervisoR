# supervisoR
An R package for the visualization of enrichment scores in hierarchical pathway networks.

To install the package, the following packages are needed:
```
BiocManager::install(c("msigdbr","igraph","ggraph","ggplot2","stringr","colorspace","ggimage",
  "magick","cowplot","magrittr","dplyr","ragg","base64enc","digest","progressr","readr","rlang",
  "shiny","visNetwork","shinyWidgets")
)
```
Then run the following to install supervisoR:
```
remotes::install_git("abmarstrand/supervisoR")
```
