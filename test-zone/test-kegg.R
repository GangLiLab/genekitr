rm(list=ls())
library(dplyr)
library(AnnoGenes)
data(geneList, package="DOSE")
id = names(geneList)[1:100]

entid = ids
ensid = transId(ids, trans_to = 'ensembl',org = 'human')%>% stringi::stri_remove_na()
symid = transId(ensid, org = 'human', trans_to = 'symbol') %>% stringi::stri_remove_na()

genKEGG(entid, org = 'hg',use_symbol = F)
genKEGG(entid, org = 'hg',use_symbol = T)

genKEGG(ensid, org = 'hg',use_symbol = F)
genKEGG(ensid, org = 'hg',use_symbol = T)
