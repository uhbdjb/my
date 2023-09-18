library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db) # mouse全基因组注释包，用于不同版本基因名的转换
library(dplyr)
library(stringr)
library(msigdbr) 

# GSEA #mdigbd dataset
msigdbr_df <- msigdbr(species = "Mus musculus", category = "C2")
# msigdbr_df$gene_symbol:"SYMBOL";$entrez_gene:"ENTREZID"
# msigdbr_list <- split(x = msigdbr_df$entrez_gene, f = msigdbr_df$gs_name)
# msigdbr_t2g <- msigdbr_df %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()
# msigdbr_t2g <- msigdbr_df %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()
msigdbr_t2g <- msigdbr_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

files <- list.files(getwd(), pattern = "*.csv") 

# 多组数据，用for循环
for (i in 1:length(files)){
  file <- files[i]
  print(file)
  dif_genes <- read.csv(file)
  # 得到差异分析结果：all_different_genes
  geneList <- dif_genes$avg_logFC
  # 如果是Ensemble ID，并且如果还带着版本号，需要去除版本号，再进行基因ID转换，得到Entrez ID
  names(geneList) <- dif_genes$gene_id
  # 最后从大到小排序，得到一个字符串
  geneList <- sort(geneList, decreasing = T) 
  
  # 检查是否有重复基因名
  genelist <- dif_genes$gene_id
  genelist <- genelist[duplicated(genelist)]
  
  gsea_result <- GSEA(geneList = geneList,
                      minGSSize = 1,
                      maxGSSize = 1000,
                      pvalueCutoff = 1,
                      TERM2GENE = msigdbr_t2g)
  # 将其中的基因名变成symbol ID
  #gsea_result <- setReadable(gsea_result, OrgDb = org.Mm.eg.db)
  # 还可以直接点击查看，只需要转成数据框
  #gsea_result_df <- as.data.frame(gsea_result)
  # 导出结果
  write.csv(gsea_result@result, file = "gsea_result.csv")
}