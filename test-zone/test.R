rm(list = ls())
options(stringsAsFactors = F)

library(stringr)
library(dplyr)
library(openxlsx)
id =c('Ticam2
Arhgap33os
Insl3
Myo15
Gal3st2b
Bloc1s1') 
id=str_split(id,"\n")[[1]]

gene_info = genInfo(id, org = 'mm')
pmd_polycomb = genPubmed(id, keywords = 'polycomb', field = 'tiab')
pmd_stemcell = genPubmed(id, keywords = 'stem cell', field = 'tiab')

wb <- createWorkbook()
wb=expo_sheet(wb, sheet_dat = gene_info, sheet_name = .nm(gene_info)) %>% 
  expo_sheet(., sheet_dat = pmd_polycomb, sheet_name = .nm(pmd_polycomb)) %>% 
  expo_sheet(., sheet_dat = pmd_stemcell, sheet_name = .nm(pmd_stemcell)) 
  

saveWorkbook(wb, "~/Downloads/test4.xlsx", overwrite = T)




