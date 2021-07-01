rm(list = ls())
options(stringsAsFactors = F)

library(stringr)
library(dplyr)
library(openxlsx)

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
# pubmed
pmd_polycomb = genPubmed(id, keywords = 'polycomb', field = 'tiab')
pmd_stemcell = genPubmed(id, keywords = 'stem cell', field = 'tiab')

# gsea
data(geneList, package="DOSE")
msigdb <- getMsigdb(org='human', category='C3',subcategory = 'TFT:GTRD')
gsea = genGSEA(genelist = geneList,geneset = msigdb,pvalueCutoff = 0.01) %>%
  DOSE::setReadable(., OrgDb = org.Hs.eg.db, keyType = 'ENTREZID') %>%
  as.data.frame()


wb <- createWorkbook()
wb=expo_sheet(wb, sheet_dat = gene_info, sheet_name = 'gene_info') %>%
  # expo_sheet(., sheet_dat = pmd_polycomb, sheet_name = 'pmd_polycomb') %>%
  expo_sheet(., sheet_dat = pmd_stemcell, sheet_name = 'pmd_stemcell') %>%
  expo_sheet(., sheet_dat = gsea, sheet_name = 'gsea')

saveWorkbook(wb, "~/Downloads/test4.xlsx", overwrite = T)




