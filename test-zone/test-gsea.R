rm(list = ls())
data(geneList, package="DOSE")
head(geneList)

source('./R/genGSEA.R')
source('./R/utils.R')
source('./R/genInfo.R')

msigdb <- getMsigdb(org='human', category='C5',subcategory='C5')


sublist=head(geneList,100)
names(sublist) = genInfo(names(sublist), org = 'hs') %>%
  dplyr::pull(symbol)

msigdb_category_data()

egmt <- genGSEA(genelist = sublist, species = 'human',
        category = 'C5',subcategory='CP:KEGG' ,
        pvalueCutoff=1,minGSSize = 5,maxGSSize = 5000)
head(egmt)

gene <- names(geneList)[abs(geneList) > 2]

egmt <- enricher(gene, TERM2GENE=msigdb)

## 来自：https://support.bioconductor.org/p/122496/
# 意思是：GSEA需要利用全部的基因集（附带logFC）进行分析【当然我们的函数支持symbol或entrez】；
# 如果只对某部分基因感兴趣，可以用ORA（如GO和KEGG去分析）
# - GSEA works on the full gene vector, testing whether genes of a gene set (here: a KEGG pathway) accumulate at the top or bottom of the full gene vector ordered by direction and magnitude of expression change.
# - If you are interested only in genes of significant expression change (here: abs(log2FC) > 2), and want to know whether certain gene sets (here: KEGG pathways) contain a disproportional number of these significant genes, you would rather carry out a over-representation analysis (ORA).
# - In clusterProfiler: use enrichKEGG for that purpose.

## 来自：https://www.biostars.org/p/9468150/
# - enrichr() performs a hypergeometric test comparing the set of "significant" genes against the "universe" (or background) genes.
# - GSEA() is a Komolgorov-Smirnov test on the whole ranked gene list, testing if some category (e.g., a KEGG pathway) is more abundant at the top of the list than expected by chance.
# - GSEA is generally considered preferable, as 1) it is considered more sensitive in general, as it can detect small, but concerted expression changes, and 2) it doesn't need the somewhat subjective step of defining which are the "significant" genes.







