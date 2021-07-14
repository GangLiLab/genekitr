rm(list = ls())
library(AnnoGenes)
# use huamn id
data(geneList, package = 'DOSE')
id = names(geneList)[1:5]
id

genInfo(id,org = 'human')

id2 = transId(id, trans_to = 'symbol',org='hs', return_dat = F)
genInfo(id2,org = 'human')

id3 = transId(id, trans_to = 'ens',org='hs', return_dat = F)
genInfo(id3,org = 'human')

# in this example, BCC7 is the alias of TP53; SXHFJG is a fake name
id4 = c("MCM10",  "CDC20",  "S100A9", "FOXM1",  "KIF23",  "MMP1",   "CDC45",  "BCC7" ,  "SXHFJG", "TP53"  )
genInfo(id4,org = 'human')

# change organism: moouse
mm_id =c('Ticam2
Arhgap33os
Insl3
Myo15
Gal3st2b
Bloc1s1
TP53')
mm_id=stringr::str_split(mm_id,"\n")[[1]]
genInfo(mm_id,org = 'mouse')


# change organism: dm fly
library(org.Dm.eg.db)
fly_id = toTable(org.Dm.egSYMBOL) %>% dplyr::pull(1) %>% sample(10)
fly_id = c(fly_id,'1') # add a human fake id
x = genInfo(fly_id,org = 'fly')


