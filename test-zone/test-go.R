rm(list=ls())
library(dplyr)
library(ggplot2)
library(genekitr)
data(geneList, package="genekitr")
id = names(geneList)[1:100]

entid = id
symid = transId(id, trans_to = 'symbol',org = 'human')
ensid = transId(id, trans_to = 'ensembl',org = 'human')

ego = genGO(entid, org = 'human',ont = 'bp',pvalueCutoff = 0.05,qvalueCutoff = 0.05 ,use_symbol = F)
ego2 = genGO(entid, org = 'human',ont = 'bp',pvalueCutoff = 0.05,qvalueCutoff = 0.05 ,use_symbol = T)

ego3 = genGO(entid, org = 'human',ont = 'cc',pvalueCutoff = 0.05,qvalueCutoff = 0.05 ,use_symbol = T)





