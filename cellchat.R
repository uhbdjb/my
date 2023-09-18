#cellchat
#首先得清空环境
merge<-merge(x=macro3,y=fib2)
merge <- NormalizeData(merge, normalization.method = "LogNormalize", scale.factor = 10000) 
DefaultAssay(object = merge) <- "RNA" 
merge <- FindVariableFeatures(merge, selection.method = "vst", nfeatures = 6000) 
merge <- ScaleData(merge)
merge <- RunPCA(merge) 
merge <- FindNeighbors(merge, dims = 1:10) 
merge <- FindClusters(merge,resolution = 0.3)
merge <- RunUMAP(merge, dims = 1:10)
DimPlot(mydata)
Cells.sub.markers <- FindAllMarkers(merge,only.pos = TRUE, win.pct =0.25, logfc.threshold = 0.25)
#if(!require(dplyr))install.packages("dplyr")
Cells.sub.markers %>% group_by(cluster) %>% top_n(n =2,wt = avg_log2FC)

top10 <- Cells.sub.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Cells.sub.markers, features = top10$gene,label=F)

VlnPlot(merge, features = c('SPP1','APOE'))
VlnPlot(merge, features = c('CSF3R'),pt.size = 0,ncol=3)+NoLegend()
FeaturePlot(merge,features = c('POSTN','RGS5'),reduction = "umap",pt.size = 0.5)
rm(list = ls())
install.packages("NMF")
library(installr)
install.Rtools()
devtools::install_github('sqjin/CellChat')
pip install umap-learn
devtools::install_github(\jokergoo/ComplexHeatmap\)
BiocManager::install("navinlabcode/copykat")#下载cellchat,可能需要在R里进行，不太好下载
1suppressMessages(if(!require(CellChat))devtools::install_github('sqjin/CellChat'))
library(CellChat)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)

#需要安装rtools才能运行
options(stringsAsFactors = FALSE) #输入数据不自动转换成因子（防止数据格式错误）
options(futrue.globlas.Maxsize=2*1024**3)
#设置硬件参数，8线程
suppressWarnings(suppressMessages(future::plan("multiprocess", workers = 8)))

#下载数据
#win+r wmic cpu get numberOfLogicalProcessors
#load(url("https://ndownloader.figshare.com/files/25950872"))#读取示例数据集
#View(data_humanSkin)
#saveRDS(data_humanSkin,'data_humanSkin.rds')

class(merge)
## [1] "list"
#示例数据：来源于人类皮肤
#作者发表文献 #https://www.nature.com/articles/s41467-021-21246-9.pdf
merge$labels <- merge@active.ident
unique(merge$labels)

data.input <- merge@assays$RNA@data # normalized data matrix
meta <- merge@meta.data # a dataframe with rownames containing cell mata data
unique(meta$group)
## [1] "LS" "NL"
cell.use <- rownames(meta)[meta$group == "Tumor"] # 按指定的变量提取细胞
data.input <- data.input[, cell.use]#取出对应细胞,也就是说，data的列名是meta的行名
meta = meta[cell.use, ]#取出对应细胞的meta信息
unique(meta$labels)#看meta中储存的细胞注释信息，稍后用它作为分组依据

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
#创建celllchat对象，group.by指定通讯间的对象，用meta中的注释作为分组依据
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
groupSize <- as.numeric(table(cellchat@idents)) #每种细胞的细胞数量

# 选择合适的物种，可选CellChatDB.human, CellChatDB.mouse
CellChatDB <- CellChatDB.human
#查看数据库的组成比例
showDatabaseCategory(CellChatDB)
# 查看数据库具体信息
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")#取出相应分类用作分析数据库
cellchat@DB <- CellChatDB.use#将数据库内容载入cellchat对象中，相当于设置好接下来要参考的数据库

cellchat <- subsetData(cellchat)#取出表达数据
cellchat <- identifyOverExpressedGenes(cellchat)#寻找高表达的基因#
cellchat <- identifyOverExpressedInteractions(cellchat)#寻找高表达的通路
cellchat <- projectData(cellchat, PPI.human)#投影到PPI，储存上一步的结果到cellchat@LR$LRsig

#计算细胞与细胞之间通信的概率
cellchat <- computeCommunProb(cellchat, raw.use = T)#默认计算方式为#type = "truncatedMean",
##去掉通讯数量很少的细胞，默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat <- filterCommunication(cellchat, min.cells = 10)

#将推断的结果提取出来
df.net <- subsetCommunication(cellchat)#将细胞通讯预测结果以数据框的形式取出
#install.packages("DT"),DT就是展示数据框的一种工具，可以调节展示的参数等
library(DT)
DT::datatable(df.net)
write.csv(df.net,'02.df.net.csv')

#计算信号通路水平上的细胞间通讯，通过汇总与每个信号通路相关的所有配体-受体相互作用的通信概率来计算信号通路水平的通信概率
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)#计算细胞-细胞聚合通信网络
groupSize <- as.numeric(table(cellchat@idents))
groupSize
#组图，1行3列
par(mfrow = c(1,3), xpd=TRUE)
#使用圆图显示任意两个细胞群之间的相互作用数量或总相互作用强度(权重)。颜色和source是一致的，圈圈的大小是每个细胞群细胞的数量
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")
#总相互作用强度
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
#以cDC2为中心的相互作用强度
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 targets.use = 'SPP1_Macro')