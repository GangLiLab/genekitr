rm(list=ls())
library(dplyr)
library(AnnoGenes)
data(geneList, package="DOSE")
ids = names(geneList)[1:100]
gkeg <- genKEGG(ids, org = 'hg',use_symbol = T)
head(gkeg)
# mapId(id = ids, from = 'entrez', to = 'symbol',org='hg', return_dat = T)
# mapId(id = ids, from = 'entrez', to = 'symbol',org='human', return_dat = T)

# if we use symbol gene as input:
test = mapId(id = names(geneList)[100:300], from = 'entrez', to = 'symbol',org='hs', return_dat = F)
head(test)
keg_raw <- genKEGG(test, org = 'hs', use_symbol = F)
keg_readable <- genKEGG(test, org = 'hs', use_symbol = T)

