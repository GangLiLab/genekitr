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

test1 = genInfo(id, org = 'mm')
test2=genPubmed(id, keywords = 'stem cell', field = 'tiab')

wb <- createWorkbook()
wb=expo_sheet(wb, sheet_name = 'genInfo',sheet_dat = test1) %>% 
  expo_sheet(., sheet_name = 'genPub',sheet_dat = test2)

saveWorkbook(wb, "~/Downloads/test3.xlsx", overwrite = T)




