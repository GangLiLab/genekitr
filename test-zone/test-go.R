rm(list=ls())
library(dplyr)
library(AnnoGenes)
data(geneList, package="DOSE")
id = names(geneList)[1:100]

ego <- genGO(id, org = 'bt',ont = 'CC',pvalueCutoff = 0.05,qvalueCutoff = 0.2)
head(ego)

