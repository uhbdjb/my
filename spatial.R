library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

setwd('~/workspace/projects/ContentShare/spatial_transcriptome')
library(SeuratData)
InstallData("stxBrain")


# Keys： spatial analysis in Seurat
## load data
## data structure 
## new function: SpatialFeaturePlot， SpatialDimPlot, LinkedDimPlot
## subset datasets or images


## Load Your Own Dataset
lymph.img = Seurat::Read10X_Image('./data/Human_Lymph_Node/spatial', image.name = 'tissue_lowres_image.png')
lymph = Seurat::Load10X_Spatial(data.dir='./data/Human_Lymph_Node',
                                filename = 'V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5',
                                assay='Spatial',
                                slice='slice1',
                                image = lymph.img
)


## 数据加载及归一化
brain <- LoadData("stxBrain", type = "anterior1")  # understand the data structure HERE
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial", slot='data') + theme(legend.position = "right")
wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
DefaultAssay(brain)

## 展示基因
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr")) #  + theme(legend.position = "right")
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 3)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2

SpatialFeaturePlot(brain, features = c('Ttr'), pt.size.factor = 5, alpha =c(0.1, 1.0))

## 降维、聚类
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3, group.by = 'seurat_clusters')
p1 + p2

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain,
                                                          idents = c(2, 1, 4, 3, 5, 8)), 
               facet.highlight = TRUE, ncol = 3)

## 交互式操作
SpatialDimPlot(brain, interactive = TRUE)
SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)
LinkedDimPlot(brain)

## 相对于交互式操作，直接subset也许更快更方便
SpatialDimPlot(brain_subset, crop = F)
brain_subset = subset(x=brain, subset=seurat_clusters==c(5, 11))
SpatialDimPlot(brain_subset, group.by = 'seurat_clusters',
               pt.size.factor = 2.0, crop = F)

## 差异基因分析(该分析作用不大，根据科学问题来做选择)
#### 1. 基于亚群的差异分析
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3],
                   alpha = c(0.1, 1), ncol = 3)

#### 2. 基于空间位置的差异分析，类似于spatialDE
# selection.method = c("markvariogram", "moransi")
brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT",
                                       features = VariableFeatures(brain)[1:50],  # 以前50个基因为例
                                       selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 3)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))
var_features = brain@assays$SCT@meta.features
var_features = var_features[VariableFeatures(brain)[1:50], ]

## 与单细胞RNA-seq数据整合分析
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
SpatialDimPlot(cortex, crop = T, label = TRUE)
allen_reference <- readRDS("./data/allen_cortex.rds")
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex,
                               reference.assay = 'SCT', query.assay = 'SCT',
                               , normalization.method = "SCT", verbose = T, reduction = 'cca')
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))
### 多切片，多样本分析
DefaultAssay(brain) <- 'SCT'
brain2 <- LoadData("stxBrain", type = "posterior1")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)
brain.merge <- merge(brain, brain2)

DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)

DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
SpatialDimPlot(brain.merge)
SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))
SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"), images = 'anterior1')


