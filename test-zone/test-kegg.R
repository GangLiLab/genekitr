rm(list=ls())
library(dplyr)
library(genekitr)
data(geneList, package="genekitr")
ids = names(geneList)[1:100]

entid = ids
ensid = transId(ids, trans_to = 'ensembl',org = 'human')
symid = transId(ensid, org = 'human', trans_to = 'symbol') %>% stringi::stri_remove_na()

genKEGG(entid, org = 'hg',use_symbol = F)
genKEGG(symid, org = 'hg',use_symbol = T)
genKEGG(ensid, org = 'hg',use_symbol = F)

