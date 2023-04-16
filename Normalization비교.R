#Normalization 방식 비교

library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)
options(future.globals.maxSize = 1e+09)

brain <- LoadData("stxBrain", type = "anterior1")
brain[["Spatial"]] <- as(brain[["Spatial"]], Class = "Assay5")

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

brain <- SCTransform(brain, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
brain <- NormalizeData(brain, verbose = FALSE, assay = "Spatial")

brain <- GroupCorrelation(brain, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
brain <- GroupCorrelation(brain, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

p1 <- GroupCorrelationPlot(brain, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- GroupCorrelationPlot(brain, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p1 + p2