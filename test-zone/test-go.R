rm(list=ls())
source('./R/genGO.R')
source('./R/utils_process.R')
source('./R/utils_analyse.R')
data(geneList, package="DOSE")
id = names(geneList)[1:100]

ego <- genGO(id, org = 'human',ont = 'CC',pvalueCutoff = 0.05,qvalueCutoff = 0.2)
head(ego)

