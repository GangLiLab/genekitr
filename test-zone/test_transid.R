library(AnnoGenes)
#--- huamn id ---#
data(geneList, package = 'DOSE')
id = names(geneList)[1:10]
transId(id, trans_to = 'symbol',org='hs', return_dat = T)
transId(id, trans_to = 'ens',org='human', return_dat = F)

# 如果选择物种不对，会提示报错
transId(id, trans_to = 'sym',org='human', return_dat = F)


#--- fly id ---#
# auto_install('org.Dm.eg.db')
library(org.Dm.eg.db)
id = toTable(org.Dm.egSYMBOL) %>% dplyr::pull(1) %>% sample(20)
id
transId(id, trans_to = 'symbol',org='fly', return_dat = F)
