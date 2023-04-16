#R, Rstudio, Rtools42 download

#Github

#Seurat Download

install.packages("githubinstall")
install.packages("devtools")
devtools::install_github("hoxo-m/githubinstall")
devtools::install_github("mangothecat/remotes")                  


remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)

install.packages("tidyverse")

devtools::install_github("thomasp85/patchwork")

install.packages("sctransform")
#remotes::install_github("satijalab/sctransform", ref="develop")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DelayedArray")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("glmGamPoi")