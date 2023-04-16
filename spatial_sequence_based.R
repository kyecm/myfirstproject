#library 불러오기
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)
options(future.globals.maxSize = 1e+09)

#stxBrain data 다운받기
InstallData("stxBrain")

#brain이라는 data frame에 자료 넣기
brain <- LoadData("stxBrain", type = "anterior1")
brain[["Spatial"]] <- as(brain[["Spatial"]], Class = "Assay5")

plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

#Normalization
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

#Gene expression visualization
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

#jpeg file로 저장하기
library(ggplot2)
plot <- SpatialFeaturePlot(brain, features = c("Ttr")) + theme(legend.text = element_text(size = 0),
                                                               legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
jpeg(filename = "C:/Users/kyecm/scRNAseq/myfirstproject/spatial_vignette_ttr.jpg", height = 700, width = 1200, quality = 50)
print(plot)
dev.off()

#spot의 크기와 점의 투명도 조절
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 2.0)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(1, 10))
p1 + p2

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,
                                                                                     5, 8)), facet.highlight = TRUE, ncol = 3)

SpatialDimPlot(brain, interactive = TRUE)

SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)

de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
