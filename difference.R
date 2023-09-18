merge$group <- merge$my.va2
unique(merge$cell_type)

table(Cells.sub$celltype)
library(plotrix)
library(dplyr)
library(ggsci)
mynames <-   table(mydata$celltype) %>% names()
myratio <-  table(mydata$celltype) %>% as.numeric()
pielabel <- paste0(mynames," (", round(myratio/sum(myratio)*100,2), "%)")


#箱线图
cellnum <- table(Cells.sub$orig.ident,Cells.sub$cell_type)
cellnum[1:5,1:5]
for (i in 1:nrow(cellnum)) {
  cellnum[i,] <- cellnum[i,]/sum(cellnum[i,])  
}
cellnum <- as.data.frame(cellnum)

library(reshape2)
colnames(cellnum) <- c('Sample','Celltype','Freq')
write.csv(cellnum,file='1.csv')
cellnum<-read.csv(file='1.csv')

cellnum$Group <- factor(cellnum$Group,levels = c('Normal','Tumor'))

unique(cellnum$Group)
library(tidyverse)

library(ggplot2)
library(ggsci)
library(ggsignif)
unique(cellnum$Group)
compaired <- list(c('Normal','Tumor'))

myplot<- list()
for (i in 1:length(unique(cellnum$Celltype))) {
  myplot[[i]] <-ggplot(data=cellnum[cellnum$Celltype==unique(cellnum$Celltype)[i],],
                       aes(x=Group,y=Freq,fill=Group))+
    geom_boxplot(width=0.2,
                 position = position_dodge(0.9))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.background=element_rect(fill='transparent', color='black'),
          axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1))+ 
    geom_point(pch=21, fill="black", color="deepskyblue")+
    scale_fill_manual(values = pal_npg("nrc")(10))+facet_wrap(~Celltype)+
    # stat_compare_means(aes(group=Group))+
    geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = 't.test')#t.test, wilcox.test
}
p.box <- Seurat::CombinePlots(myplot,ncol = 3,legend='right')
p.box

