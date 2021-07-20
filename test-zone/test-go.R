rm(list=ls())
library(dplyr)
library(ggplot2)
library(AnnoGenes)
data(geneList, package="DOSE")
id = names(geneList)[1:20]

entid = id
symid = transId(id, trans_to = 'symbol',org = 'human')
ensid = transId(id, trans_to = 'ensembl',org = 'human')

ego = genGO(entid, org = 'human',ont = 'bp',pvalueCutoff = 0.01,qvalueCutoff = 0.01 ,use_symbol = F)
genGO(entid, org = 'human',ont = 'bp',pvalueCutoff = 0.01,qvalueCutoff = 0.01 ,use_symbol = T)

genGO(symid, org = 'human',ont = 'bp',pvalueCutoff = 0.01,qvalueCutoff = 0.01)

genGO(ensid, org = 'human',ont = 'bp',pvalueCutoff = 0.01,qvalueCutoff = 0.01, ,use_symbol = F)
genGO(ensid, org = 'human',ont = 'bp',pvalueCutoff = 0.01,qvalueCutoff = 0.01, ,use_symbol = T)



