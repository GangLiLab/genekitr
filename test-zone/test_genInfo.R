rm(list = ls())
library(AnnoGenes)
#--- huamn id ---#
data(geneList, package = 'DOSE')
id = names(geneList)[1:5]
id

genInfo(id,org = 'human')
