library(Seurat)
library(dplyr)
setwd("xx")


ctrl2 <- Seurat::Read10X(data.dir = "./Control2")
ad1 <- Seurat::Read10X(data.dir = "./AD1")
seu_Control <- Seurat::CreateSeuratObject(counts = ctrl2,
                                  project = "p2",
                                  min.cells = 3,
                                  min.features = 100)
seu_AD1 <- Seurat::CreateSeuratObject(counts = ad1,
                                  project = "p2",
                                  min.cells = 3,
                                  min.features = 100)
exp <- rep("ctrl",dim(seu_Control@meta.data)[1])
seu_Control@meta.data <-cbind(seu_Control@meta.data,exp)
exp <- rep("ad",dim(seu_AD1@meta.data)[1])
seu_AD1@meta.data <-cbind(seu_AD1@meta.data,exp)

seu.list <- list(seu_Control,seu_AD1)
seu.list <- lapply(X = seu.list, FUN = function(x){
                    x <- NormalizeData(x)
                    x<- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})
features <- SelectIntegrationFeatures(object.list = seu.list)
anchors <- FindIntegrationAnchors(object.list = seu.list, anchor.features = features)  #warning, same cell name
seu <- IntegrateData(anchorset = anchors)
DefaultAssay(seu) <- "integrated"
# seu <- merge(seu_Control,seu_AD1)
# table(seu@meta.data$exp)
# seu <- Seurat::PercentageFeatureSet(seu, 
#                                     pattern = "^MT-", 
#                                     col.name = "percent.mito")
# seu <- subset(seu, subset = nFeature_RNA > 200 & 
#                 nFeature_RNA < 5000 &
#                 percent.mito < 10)
# seu <- Seurat::NormalizeData(seu,
#                      normalization.method = "LogNormalize",
#                      scale.factor = 10000)
# seu <- Seurat::FindVariableFeatures(seu,
#                             selection.method = "vst",
#                             nfeatures = 2000)
seu <- Seurat::ScaleData(seu,
                 features = rownames(seu))
# top <- head(Seurat::VariableFeatures(seu), 5)
# vf_plot <- Seurat::VariableFeaturePlot(seu)
# Seurat::LabelPoints(plot = vf_plot,
#             points = top10, repel = TRUE)

seu <- Seurat::RunPCA(seu, npcs = 25, verbose = FALSE)
# Seurat::DimPlot(seu, reduction = "pca")
# Seurat::ElbowPlot(seu, ndims = 40)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:25)
Seurat::DimPlot(seu, reduction = "umap")

seu <- Seurat::FindNeighbors(seu, dims = 1:25)
seu <- Seurat::FindClusters(seu, resolution = seq(0.1, 0.8, by=0.1))
head(seu@meta.data)
library(clustree)
clustree::clustree(seu@meta.data[,grep("integrated_snn_res", colnames(seu@meta.data))],
                   prefix = "integrated_snn_res.")   #res.
Seurat::DimPlot(seu, group.by = "integrated_snn_res.0.7")
seu <- Seurat::SetIdent(seu, value = seu$integrated_snn_res.0.7)
DefaultAssay(seu) <- "RNA"

# library(celldex)
# library(SingleR)
# ref <- celldex::BlueprintEncodeData()   ##bfc err epuf, code sib NAse
# class(ref)
# table(ref$label.main)
# ref_reduce <- ref[,ref$label.main %in% c("B-cells","CD4+ T-cells","CD8+ T-cells","Monocytes","NK cells","Epithelial cells","HSC","Eosinophils","DC")]
# seu_SingleR <- SingleR(test=GetAssayData(seu, assay = "RNA", slot = "data"),
#                            ref=ref_reduce, labels=ref_reduce$label.fine)
# head(seu_SingleR)
# SingleR::plotScoreHeatmap(seu_SingleR)
# seu$SingleR_annot <- seu_SingleR$labels
# dittoSeq::dittoDimPlot(seu, "SingleR_annot", size = 0.7)
# #all Epithelial

de_genes <- Seurat::FindAllMarkers(seu,  min.pct = 0.25,
                                   only.pos = TRUE)
tops <- subset(de_genes, de_genes$p_val_adj<0.05)
write.csv(tops,
          "de_genes_FindAllMarkers.csv",
          row.names = F, quote = F)
top_specific_markers <- de_genes %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)
write.csv(top_specific_markers,
          "top_specific_markers.csv",
          row.names = F, quote = F)
# tops <- Seurat::FindMarkers(seu,
#                            group.by=seu@meta.data$exp,
#                            ident.1="ad",ident.2="ctrl")
annotations <- read.csv("annotation.csv")
# install.packages('Rtools')
# install.packages('metap')
# "metap" %in% installed.packages().  #true
# BiocManager::install('multtest')
# library(metap)
table(seu@meta.data$integrated_snn_res.0.7)
table(seu@meta.data$exp)
lapply(seu@meta.data, function(x) as.data.frame(table(x, seu@meta.data$exp)))
#https://stackoverflow.com/questions/60824490/how-to-run-the-table-function-in-r-across-multiple-variables-and-compile-the-res
saveRDS(seu, "seu-int.rds")
readRDS("seu-int.rds")


# nk.markers <- Seurat::FindConservedMarkers(seu, ident.1 = 0, grouping.var = "exp", verbose = FALSE)    #Seurat::

cluster0_conserved_markers <- FindConservedMarkers(seu,
                              ident.1 = 0,
                     	      grouping.var = "exp",
                              only.pos = TRUE,
		              logfc.threshold = 0.25)
library("tidyr")
library(tibble)
cluster0_ann_markers <- cluster0_conserved_markers %>%
                rownames_to_column(var="gene") %>%
                left_join(y = unique(annotations[, c("gene_name", "description")]),
                          by = c("gene" = "gene_name"))
View(cluster0_ann_markers)

get_conserved <- function(cluster){
  FindConservedMarkers(seu,
                       ident.1 = cluster,
                       grouping.var = "exp",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
               by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
  }
conserved_markers <- purrr::map_dfr(c(0,1,2,3,4,5,6,7,8), get_conserved)
top10 <- conserved_markers %>% 
  mutate(avg_fc = (ctrl_avg_log2FC + ad_avg_log2FC)/2 ) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, wt = avg_fc)
View(top10)
write.csv(top10,
          "top10-conserved.csv",
          row.names = F, quote = F)
FeaturePlot(object = seu,features = c("HIST1H1B", "HIST1H1A", "HIST1H1D", "HIST2H2AC"),
                         sort.cell = TRUE,
                         min.cutoff = 'q10', 
                         label = TRUE,
			                   repel = TRUE)
VlnPlot(seu, features = c("HIST1H1B", "HIST1H1A", "HIST1H1D", "HIST2H2AC"))
seurat_integrated <- RenameIdents(object = seu, 
                               "0" = "CD71+ EarlyErythroid",
                               "1" = "B lymphoblasts",
                               "2" = "CD8+ Proliferating T",
                               "3" = "Myofibroblast",
                               "4" = "CD4+ Proliferating T",
                               "5" = "Stromal cells",
                               "6" = "CD14+ Monocytes",
                               "7" = "Erythroblasts/Vascular endothelial cells",
                               "8" = "Natural Killer T (NKT) cell")
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)
saveRDS(seu, "seu.rds")
seu <- readRDS("seu.rds")

library(limma)
library(edgeR)

#xxx, below

ad-vs-control-0 <- Seurat::FindMarkers(seu0,)
ad-vs-control-0 <- subset(ad-vs-control-0, ad-vs-control-0$p_val_adj<0.05)
# merge_limma_FindMarkers <- merge(cluster6_vs_cluster0, limma_de, by="row.names", all.x=T)
# par(mar=c(4,4,4,4))
# plot(merge_limma_FindMarkers$avg_log2FC, merge_limma_FindMarkers$logFC, xlab="log2FC Wilcoxon", ylab="log2FC limma", pch=15, cex=0.5)
# abline(a=0, b=1, col="red")




seu0 <- subset(seu, idents = 8)
sce_limma <- as.SingleCellExperiment(seu0)
y <- DGEList(counts(sce_limma), samples=colData(sce_limma))   #edgeR
keep <- filterByExpr(y, group=y$samples$SingleR_annot)
y <- y[keep,]
summary(keep)
#    Mode   FALSE    TRUE 
# logical   18200     168 
y$samples$Clustering <- y$samples$ident  #res0.7
y$samples$exp <- y$samples$exp
design <- model.matrix(~0 + exp, data = y$samples)
design
contrast.mat <- limma::makeContrasts(expctrl, levels = design)
dge <- edgeR::calcNormFactors(y) 
vm <- limma::voom(dge, design = design, plot = TRUE)
fit <- limma::lmFit(vm, design = design)
fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
fit.contrasts <- limma::eBayes(fit.contrasts)
limma::topTable(fit.contrasts, number = 10, sort.by = "P")
limma_de <- limma::topTable(fit.contrasts, number = Inf, sort.by = "P")
length(which(limma_de$adj.P.Val<0.05))  #only tails
limma_de0 <- subset(limma_de, limma_de$adj.P.Val<0.05)
write.csv(limma_de0,
          "limma_de8",
          row.names = T, quote = F)




#check 2 annotations:`singleR` and `CHETAH` 
seu$chetah <- parotid10x.sce$celltype_CHETAH
#head(seu@meta.data)
#https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1795-z
#dittoSeq::dittoBarPlot(seu, var = "chetah", group.by = "orig.ident")
#https://r-coder.com/barplot-r/, dont need stack, can plot line, mosaic
#plot(is.factor(seu$SingleR_annot), is.factor(seu$chetah))   #mieux gsub
an_s <- sub(".*B.*cell.*", "1", seu$SingleR_annot)
an_s <- sub(".*T.*cell.*", "2", an_s)
#an_s <- sub("^T&^B", "N", an_s)
an_s <- sub("...*", "3", an_s)
an_t <- sub(".*B.*cell", "1", seu$chetah)
an_t <- sub(".*T cell.*", "2", an_t)
an_t <- sub("...*", "3", an_t)  #cuz DC

#check doublets
table(seu$Doublets == "Singlet" & seu$Supervised_Singlets == "Keep")  #cant &&
Seurat::DimPlot(seu, reduction = "umap", group.by = c("Doublets","Supervised_Singlets"))
head(seu@meta.data)
seu$sd <- ifelse((seu$Doublets == "Singlet" & seu$Supervised_Singlets == "Keep") == 1,"1","2")
Seurat::DimPlot(seu, reduction = "umap", group.by = "sd")
seut <- subset(seu, subset = sd == "1")


#cross-talk
memb <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
memb[5, ] <- mat[5, ]
netVisual_circle(memb, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Memory B-B") 
netVisual_bubble(cellchatParotid10x, sources.use = c(1:4), targets.use = 5, remove.isolate = FALSE)
cellchat <- cellchatParotid10x
cellchat@netP$pathways      #[1] "MIF"  "CD70"
pathways.show <- c("MIF")  #Cannot find  CD99
pairLR.CD70 <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CD70[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
plotGeneExpression(cellchat, signaling = "MIF")


#trajectory
library(slingshot)
rd <- reducedDim(sce, "PCA")[,1:2]
cl <- sce$seurat_clusters
names(cl) <- colnames(sce)
lin2 <- getLineages(rd, cl, start.clus= '1', end.clus = '3')  #3 red
Farben <- brewer.pal(9,"Set1")[cl]
plot(rd, col = Farben, asp = 1, pch = 16)
lines(SlingshotDataSet(lin2), lwd = 3, col = 'black', show.constraints = TRUE)
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')





