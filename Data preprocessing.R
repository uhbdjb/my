if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("celldex")
if(!require(patchwork))install.packages("patchwork")
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
nBiocManager::install("Seurat")
BiocManager::install("celldex")
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
search()
detach(celldex)
rm(list = ls())
setwd('E:/colon/GSE/my')
getwd()

library(data.table)
s = Sys.time()
new_counts <- as.matrix(fread("GSE149614_HCC.scRNAseq.S71915.count.txt"),rownames=1,header = T,stringsAsFactors = F)
e = Sys.time()
print(e-s)

library(data.table)
theUrl <- "http://www.jaredlander.com/data/TomatoFirst.csv"
new_counts <- fread(input='GSE206785_scgex.txt.gz',sep=',', header=TRUE)


## read matrix
new_counts <- read.table( "GSE149614_HCC.scRNAseq.S71915.count.txt", sep="\t", header=T, row.names=1)
?file.path
dim(matrix_data)
matrix_data[1:4,1:4]

6seurat_obj <- CreateSeuratObject(counts = matrix_data)
new_counts <- Read10X(data.dir = "cancer/") 
#读取数据
#mydata <- CreateSeuratObject(counts=mydata.data,project='pmbc3k',min.cells=3,min.features=200)
new_counts <- read.csv(file="GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv",header = TRUE, sep = ",",quote="\"", dec=".")
new_counts <- read.table(file="GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt",header = FALSE, sep = "",quote="\"", dec=".")
?read.table
new_counts<-new_counts[!duplicated(new_counts$Index), ] 
row.names(new_counts)=new_counts[,1]
new_counts=new_counts[,-1]
#回收垃圾
gc()
memory.size(T)
memory.size(F)
mydata = mydata %>% select(starts_with(c('LUNG_N01','LUNG_N06','LUNG_N08','LUNG_N09','LUNG_N20','LUNG_N28','LUNG_N30','LUNG_N31',
                                         'LUNG_T06','LUNG_T08','LUNG_T09','LUNG_T18','LUNG_T20','LUNG_T28','LUNG_T30','LUNG_T31')))
new_counts = new_counts %>% select(starts_with(c('HCC01','HCC02','HCC03','HCC04',
                                                 'HCC05','HCC08T','HCC08N')))

head(new_counts)
mydata <- CreateSeuratObject(counts = mydata,min.cells=3,project="mydata_scRNAseq")
mydata
pbmc <- CreateSeuratObject(counts = new_counts,
                           project = "pbmc3k",
                           min.cells = 3,
                           min.features = 200)
#ncol(mydata)
a <-as.data.frame(mydata[["RNA"]]@counts)
#write.table(a,"mycount.txt",sep='\t')
mydata[["percent.MT"]] <- PercentageFeatureSet(mydata,pattern="^MT-")
#人源需换成MT
head(mydata@meta.data,5)

VlnPlot(mydata,features=c("nFeature_RNA","nCount_RNA","percent.MT",ncol=3))

#if(!require(patchwork))install.packages("patchwork")
plot1<- FeatureScatter(mydata,feature1 = "nCount_RNA",feature2="percent.MT")
plot2<- FeatureScatter(mydata,feature1 = "nCount_RNA",feature2="nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))

mydata <- subset(mydata,subset = nFeature_RNA >200 & nFeature_RNA<6000 & percent.MT < 20)
#ncol(as.data.frame(mydata[["RNA"]]@counts))
mydata
mydata <- NormalizeData(mydata,normalization.method = 'LogNormalize',scale.factor = 10000)

merge <- FindVariableFeatures(merge,selection.method = "vst", nfeatures = 6000)
#for PCA DoHeatmap
top10 <-head(VariableFeatures(mydata),10)
plot1 <- VariableFeaturePlot(mydata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

all.genes <- rownames(mydata)
mydata <- ScaleData(mydata,features = all.genes)
#mydata<- ScaleData(mydata) #only variablefeatyres


mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata))
print(mydata[["pca"]],dims= 1:5, nFeatures=5)
VizDimLoadings(mydata,dims=1:2,reduction = "pca")
DimPlot(mydata, reduction = "pca")
DimHeatmap(mydata,dims=1,cells=500,balanced = TRUE)
DimHeatmap(mydata,dims=1:20,cells=500,balanced = TRUE)


mydata <- JackStraw(mydata,num.replicate = 100)
mydata <- ScoreJackStraw(mydata, dims = 1:20)
JackStrawPlot(mydata, dims = 1:20)
ElbowPlot(mydata)

mydata <- FindNeighbors(mydata, dims =1:15)
mydata <- FindClusters(mydata,resolution = 1.1)#resolution在0.4-1.2之间，越大细胞分群越多

library(clustree)
mydata <- FindClusters(
  object = mydata,
  resolution = c(seq(.1,1.6,.2))
)
clustree(mydata@meta.data, prefix = "integrated_snn_res.")

mydata.markers <- FindAllMarkers(mydata,only.pos = TRUE, win.pct =0.25, logfc.threshold = 1)
write.csv(mydata.markers,file='marker.csv')
0?#if(!require(dplyr))install.packages("dplyr")
  a<-mydata.markers %>% group_by(cluster) %>% top_n(n =2,wt = avg_log2FC)


top5 <- mydata.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(mydata, features = top5$gene,label=F,group.bar.height = 0.03, group.colors = c( "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7"))+ scale_fill_gradientn(colors = c( "navy", "white", "firebrick3"))

#多样本整合
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
library("Matrix")
rm(list = ls())
setwd("E:/colon")
getwd()

ns <- Read10X(data.dir = "ns/") 
dp <- Read10X(data.dir = "dp/")
dp <- CreateSeuratObject(counts = dp,min.cells=3,min.features=200,project="merge_scRNAseq")
ns <- CreateSeuratObject(counts = ns,min.cells=3,min.features=200,project="merge_scRNAseq")
## edit metadata
# dataset name；在seurat对象中加入project名
mydata1$group <- "a"
mydata2$group <- "b"
merge<-merge(x=mydata1,y=mydata2)
merhe <- merge(mydata1, y = c(mydata2, mydata3), add.cell.ids = c("3K", "4K", "8K"), 
               project = "PBMC15K")



# split the dataset into a list of two seurat objects (stim and CTRL)
merge.list <- SplitObject(merge, split.by = "group")

# normalize and identify variable features for each dataset independently
merge.list <- lapply(X = merge.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 6000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = merge.list)
mergegroup <- FindIntegrationAnchors(object.list = merge.list, anchor.features = features)

# this command creates an 'integrated' data assay
merge <- IntegrateData(anchorset = mergegroup)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(merge) <- "integrated"

# Run the standard workflow for visualization and clustering
merge <- ScaleData(merge, verbose = FALSE)
merge <- RunPCA(merge, npcs = 30, verbose = FALSE)
merge <- RunUMAP(merge, reduction = "pca", dims = 1:30)
merge <- FindNeighbors(merge, reduction = "pca", dims = 1:30)
merge <- FindClusters(merge, resolution = 0.5)




# Visualization
p1 <- DimPlot(Cells.sub, reduction = "umap", group.by = "group")
p2 <- DimPlot(Cells.sub, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
p1

DimPlot(Cells.sub, reduction = "umap", group.by = "sample")


