setwd('E:/colon2/spatial')
library(devtools)

DefaultAssay(lymph) <- "RNA"
cd_features <- list(c(
  'CXCL5',
  'INHBA',
  'MSC',
  'PLIN2',
  'SDC4',
  'LGALS1',
  'FN1',
  'FTL',
  'NMB',
  'S100A11',
  'IFI6',
  'C15orf48',
  'MMP9',
  'S100A10',
  'ADM',
  'CTSB',
  'CCL3',
  'MT1E',
  'ANXA2',
  'MT2A',
  'FLNA',
  'TGM2',
  'MMP14',
  'CXCL8',
  'MMP7',
  'ISG15',
  'TUBA1C',
  'SERPINE1',
  'COL1A1',
  'CTHRC1',
  'PDPN'
))
  


Inscore <- AddModuleScore(lymph,
                          features = cd_features,
                          ctrl = 100,
                          name = "CD_Features")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[8] <- 'Angiogenesis_Score' 


library(ggplot2)
mydata<- FetchData(Inscore,vars = c("UMAP_1","UMAP_2","PDPN+CAF_Score"))
a <- ggplot(coor,aes(x = UMAP_1,y =UMAP_2,colour = PDPN+CAF_Score))+
geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
colours = c('#443C85',"#2A8D86",'#8AD34E','#F9EA22'))