#生存分析
setwd("E:/colon2/survival")
#安装加载R包
#install.packages("survival")
#install.packages("forestplot")
library(survival)
library(forestplot)
library(tidyverse)
#下载生存信息
#xena官网：https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Liver%20Cancer%20(COAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#读取生存信息tsv文件
surv = read.table(file = 'TCGA-COAD.survival.tsv', sep = '\t', header = TRUE) 
#整理生存信息数据
surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]
colnames(surv) <- c('futime','fustat')
colnames(s2) <- c('fustat','futime')
#读取表达数据
expr <- read.table("COAD_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]
#表达数据整理完毕

#整合
diff<-read.csv(file='2.csv')
rownames(diff)=diff[,1]
deg_expr <- expr[rownames(diff),] %>% t() %>% as.data.frame()

deg_expr1<-scale(deg_expr)
deg_expr1<-as.data.frame(deg_expr1)
surv.expr <- cbind(s2,deg)
write.csv(surv.expr,file='d.csv')

diff<-read.csv(file='d.csv',row.names = 1)
diff$futime=diff$futime/365 

rt<-diff
####根据基因高低组做生存分析####
# 1. Determine the optimal cutpoint of variables
library(survminer)
res.cut <- surv_cutpoint(rt, #数据集
                         time = "futime", #生存状态
                         event = "fustat", #生存时间
                         variables = c('d') #需要计算的数据列名
)

summary(res.cut) #查看数据最佳截断点及统计量
# 2. Plot cutpoint for DEPDC1
# 以DEPDC1为例
plot(res.cut, "d", palette = "npg")
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(futime, fustat) ~d, data = res.cat)




ggsurvplot(fit,
           data = rt,
           pval = TRUE,
           risk.table = F, # 显示风险表
           risk.table.col = "strata",
           palette = c("#ED0000E5","#00468BE5"), # 配色采用jco
           legend.labs = c('Myeloid_High(n=84)', 'Myeloid_Low(n=385)'), # 图例
           size = 1,
           xlim = c(0,10), # x轴长度，一般为0-10年
           break.time.by = 2, # x轴步长为20个月
           legend.title = "",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Overall Survival", # 修改y轴标签
           xlab = "Time (years)", # 修改x轴标签
           ncensor.plot = F, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)

