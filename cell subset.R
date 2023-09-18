library("multtest")
library("Seurat")
library("dplyr")
library("patchwork")
library("R.utils")
library("tidyverse")
library("ggplot2")
library("SingleR")
library("MatrixGenerics")
library("celldex")

mydata<-mydata(x=Cells.sub,y=my)
mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 10000) 
DefaultAssay(object = mydata) <- "RNA" 
mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 6000) 
mydata <- ScaleData(merge1)
mydata <- RunPCA(merge1) 
mydata <- FindNeighbors(mydata, dims = 1:15) 
mydata <- FindClusters(mydata,resolution = 0.5)
mydata <- RunUMAP(mydata, dims = 1:20)

library(clustree)
merge <- FindClusters(
  object = colon2,
  resolution = c(seq(.1,1.6,.2))
)
clustree(merge@meta.data, prefix = "integrated_snn_res.")

new.cluster.ids <-c("T cells","T cells","B cells","B cells","Epithelial cells",
                    "T cells","T cells","T cells","T cells","Myeloid cells",
                    "B cells","Fibroblasts","Epithelial cells","Epithelial cells","T cells",
                    "Endothelial cells","Myeloid cells")
names(new.cluster.ids) <- levels(mydata)
mydata <- RenameIdents(mydata,new.cluster.ids)
DimPlot(mydata,reduction = 'umap')

