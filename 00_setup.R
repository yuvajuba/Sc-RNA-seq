# ================================================
#   Thesis project Setup - scRNAseq Analysis
# ================================================

## 1. Load or Install Packages
load_or_install <- function(pkgs){
  for(pkg in pkgs){
    if (!require(pkg, character.only = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
      library(pkg, character.only = TRUE)
    }
  }
}

# -- List of packages for this project
List_packages <- c("dplyr", "tibble", "stringr", "ggplot2", "tidyr", "ggsci", "shiny",
                   "ggrepel", "Seurat", "writexl", "knitr", "clustree", "conflicted",
                   "clusterProfiler", "org.Hs.eg.db", "ComplexHeatmap", "bslib")

# -- Load all required packages
load_or_install(List_packages)

# -- Set preferred packages in conflicts
prefer_dplyr <- c("filter","select","arrange","count","setdiff","intersect", "mutate")
lapply(prefer_dplyr, function(f) conflict_prefer(f, "dplyr"))

conflict_prefer("str_detect", "stringr")


## 2. Load custom project functions
function_files <- list.files("d:/Projects/JFPeyron/00_functions/", 
                             full.names = TRUE, 
                             pattern = "\\.R$")

# function_files <- list.files("~/Bureau/Projects/Thesis/00_functions", 
#                              full.names = TRUE, 
#                              pattern = "\\.R$")

sapply(function_files, source)


## 3. Define a colour palette
ColorPalettes <- list(Personal = c("#7744ff","#882200","#228800","#ffcc22","#002288","#aa0077",
                                   "#77aa00","#ff0000","#00aaff","#9900ff","#00ff00","#0000ff",
                                   "#111111","#ff7700","#00ff66","#ff22ff","#33ffff","#ffff22",
                                   "#7722cc","#003333","#00dddd","#dd33cc","#661100","#ddaa00",
                                   "#3333cc","#0077ff","#22ffaa","#44ff11","#ff00cc","#006600"),
                      ggsci_jco = ggsci::pal_jco()(10),
                      ggsci_bmj = ggsci::pal_bmj()(9),
                      ggsci_futurama = ggsci::pal_futurama()(12),
                      ggsci_cosmic = ggsci::pal_cosmic(palette = "hallmarks_light")(10))

MyPalette <- c(ColorPalettes$Personal, 
               ColorPalettes$ggsci_bmj, 
               ColorPalettes$ggsci_cosmic, 
               ColorPalettes$ggsci_jco,
               ColorPalettes$ggsci_futurama)

rm(List_packages, load_or_install, function_files, prefer_dplyr)










