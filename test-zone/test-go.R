rm(list=ls())
library(dplyr)
library(AnnoGenes)
data(geneList, package="DOSE")
id = names(geneList)[1:100]

ego <- genGO(id, org = 'human',ont = 'mf',pvalueCutoff = 0.05,qvalueCutoff = 0.1 ,readable = T)
head(ego)

plotEnrichDot(ego,xlab_type =  'FoldEnrich', legend_by = 'qvalue',show_item = 10)

