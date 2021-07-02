rm(list=ls())
library(dplyr)
library(AnnoGenes)
data(geneList, package="DOSE")
ids = names(geneList)[1:10]

mapId(id = ids, from = 'entrez', to = 'symbol',org='human', return_dat = F)

gkeg <- genKEGG(ids, org = 'human')
head(gkeg)
