library("airway")
data(airway)
airway
as.data.frame(colData(airway))
summary(colSums(assay(airway))/1e6)
metadata(rowRanges(airway))

exprSet=assay(airway)
group_list=colData(airway)[,3]
colnames(exprSet)
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

length(unique(rownames(exprSet)))
id = transId(rownames(exprSet),'symbol','hs',simple = F)

id2 = clusterProfiler::bitr(rownames(exprSet),'ENSEMBL','SYMBOL','org.Hs.eg.db')

table(is.na(id))


x = genInfo(rownames(exprSet),'hs')


