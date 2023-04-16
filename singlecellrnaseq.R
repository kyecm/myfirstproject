#https://portrai.io/untitled/

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
#C:\Users\kyecm\scRNAseq\myfirstproject\pbmc>tar -xf pbmc3k_filtered_gene_bc_matrices.tar.gz -> tar 압축 해제
pbmc.data <- Read10X(data.dir = "C:/Users/kyecm/scRNAseq/myfirstproject/pbmc/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

dense.size <- object.size(as.matrix(pbmc.data))
dense.size
