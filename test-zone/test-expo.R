rm(list = ls())
options(stringsAsFactors = F)

library(stringr)
library(dplyr)

#--- test id ---#
id =c('Ticam2
Arhgap33os
Insl3
Myo15
Gal3st2b
Bloc1s1')
id=str_split(id,"\n")[[1]]

# gene info
gene_info = genInfo(id, org = 'mm')

# # pubmed
# pmd_polycomb = genPubmed(id, keywords = 'polycomb', field = 'tiab')
# pmd_stemcell = genPubmed(id, keywords = 'stem cell', field = 'tiab')

# gsea
data(geneList, package="genekitr")
gse <- genGSEA(genelist = geneList, org = "human",
               category = "H",use_symbol = TRUE)

expoSheet(
  dat_list = list(gene_info, gse), name_list = list("gene_info", "GSEA"),
  filename = "test.xlsx", dir = '~/Downloads/')




