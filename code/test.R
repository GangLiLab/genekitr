library(stringr)
library(dplyr)
id =c('Ticam2
Arhgap33os
Insl3
Myo15
Gal3st2b
Bloc1s1') 
id=str_split(id,"\n")[[1]]

genInfo(id, org = 'mm', html_result = F, destdir = '~/Downloads/')

res3=genPubmed(id, keywords = 'stem cell', field = 'tiab')






