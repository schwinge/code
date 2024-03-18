
library(Seurat) # load the correct Seurat version
library(stardust)
library(Rfast2)
library(hdf5r)
HBC1 <- readRDS("HBC1_subset.RDS")
HBC1
head(HBC1@assays$Spatial@counts)
dim(HBC1)
p1 <- SpatialPlot(HBC1, pt.size.factor = 0)
p2 <- SpatialPlot(HBC1)   #jupo
HBC1 <- SCTransform(HBC1, assay = "Spatial", verbose = FALSE)
HBC1 <- RunPCA(HBC1, assay = "SCT", verbose = FALSE)
HBC1 <- FindNeighbors(HBC1, reduction = "pca", dims = 1:30)
HBC1 <- FindClusters(HBC1, verbose = FALSE)
HBC1 <- RunUMAP(HBC1, reduction = "pca", dims = 1:30)
p1 <- DimPlot(HBC1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(HBC1, label = TRUE, label.size = 3)
p1 + p2



HBC13 <- subset(HBC1, idents = 3)
p1 <- SpatialDimPlot(HBC13, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(HBC13, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2
