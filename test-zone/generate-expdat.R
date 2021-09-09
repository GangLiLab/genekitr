rm(list = ls())
library("airway")
data(airway)
exprSet=assay(airway)
group_list=colData(airway)[,3]
colnames(exprSet)
if(F){
  pheatmap::pheatmap(cor(exprSet))
  group_list
  tmp=data.frame(g=group_list)
  rownames(tmp)=colnames(exprSet)
  # 组内的样本的相似性应该是要高于组间的！
  pheatmap::pheatmap(cor(exprSet),annotation_col = tmp)
  dim(exprSet)
  exprSet=exprSet[apply(exprSet,1, function(x) sum(x>1) > 5),]
  dim(exprSet)
  exprSet=log(edgeR::cpm(exprSet)+1)
  dim(exprSet)
  exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:500]),]
  dim(exprSet)
  M=cor(log2(exprSet+1))
  tmp=data.frame(g=group_list)
  rownames(tmp)=colnames(M)
  pheatmap::pheatmap(M,annotation_col = tmp)
  pheatmap::pheatmap(scale(cor(log2(exprSet+1))))
}
library(genekitr)
length(unique(rownames(exprSet)))
id = transId(rownames(exprSet),'entrez','hs',unique = F) %>%  na.omit()

des = function(expr,group_list){
  colData <- data.frame(row.names =colnames(expr),
                        condition=group_list)
  dds <- DESeqDataSetFromMatrix(
    countData = expr,
    colData = colData,
    design = ~ condition)
  #参考因子应该是对照组 dds$condition <- relevel(dds$condition, ref = "untrt")

  dds <- DESeq(dds)
  # 两两比较
  res <- results(dds, contrast = c("condition",rev(levels(group_list))))
  resOrdered <- res[order(res$pvalue),] # 按照P值排序
  DEG <- as.data.frame(resOrdered)
  # 去除NA值【为何存在pvalue为NA？】
  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
  DEG <- na.omit(DEG)


  DEG$gene=rownames(DEG)
  return(list(DEG,dds))
}
library(DESeq2)
res  = des(exprSet,group_list)
DEG = res[[1]]

DEG$id = transId(rownames(DEG),'entrez','hs')
DEG  = DEG[!is.na(DEG$id),]
geneList = DEG$log2FoldChange
names(geneList) = DEG$id
geneList = sort(geneList,decreasing = T)
head(geneList)
# save(geneList,file = 'data/geneList.rda')


entid = names(geneList)[1:100]
g1 = genGO(entid, org = 'human',ont = 'bp',pvalueCutoff = 0.05,qvalueCutoff = 0.05 ,use_symbol = T)
g2 = genGO(entid, org = 'human',ont = 'all',pvalueCutoff = 0.05,qvalueCutoff = 0.05 ,use_symbol = F)

g3 = genGSEA(genelist = geneList,org = 'human', category='C3',
        subcategory = 'TFT:GTRD',use_symbol = FALSE)
g32 = genGSEA(genelist = geneList,org = 'human', category='C3',
              subcategory = 'TFT:GTRD',use_symbol = T)

g4 = genKEGG(entid, org = 'human',use_symbol = T)
g5 = genKEGG(entid, org = 'human',use_symbol = F)


