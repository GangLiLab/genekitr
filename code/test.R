library(stringr)

up_gene =c('Ticam2
Arhgap33os
Insl3
Myo15
Gal3st2b
Aoc3
Krt7
Gbp4
Gm2237
Fhit
Zbp1
Cyp2c23
Bloc1s1')

up_gene=str_split(up_gene,"\n")[[1]]

down_gn = c('Col4a1
Col4a2
Smad5
Gm8300
Txndc15
Sox7')

dn_gene=str_split(down_gn,"\n")[[1]]

new_info=gene_info %>% 
  filter(symbols %in% c(up_gene,dn_gene)) %>% 
  mutate(type = case_when(symbols %in% up_gene ~ 'up',
                          symbols %in% dn_gene ~ 'down'))

file='~/Downloads/gene_info.html'
y <- DT::datatable(new_info,escape = F,rownames=F)
DT::saveWidget(y,file)






