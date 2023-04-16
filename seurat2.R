library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)
options(future.globals.maxSize = 1e+09)

InstallData("stxBrain")

brain <- LoadData("stxBrain", type = "anterior1")
brain[["Spatial"]] <- as(brain[["Spatial"]], Class = "Assay5")

plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
