rm(list = ls())
library(AnnoGenes)
#--- huamn id ---#
data(geneList, package = 'DOSE')
id = names(geneList)[1:5]
id
transId(id, trans_to = 'symbol',org='hs', return_dat = T)

fake_id = c(id,'23215326','1','2','344263475','45')
res = transId(fake_id, trans_to = 'sym',org='human', return_dat = T)
res
only_id = transId(fake_id, trans_to = 'sym',org='human', return_dat = F)
only_id
identical(fake_id, res$entrezid)

transId(na.omit(res$symbol), trans_to = 'ens',org='hs', return_dat = T)
transId(na.omit(res$symbol), trans_to = 'entrez',org='hs', return_dat = T)

# 如果选择物种不对，会提示报错
transId(id, trans_to = 'sym',org='mouse', return_dat = F)

# transId vs bitr
fake_id
transId(fake_id, trans_to = 'sym',org='human', return_dat = T)
clusterProfiler::bitr(fake_id, fromType = 'ENTREZID',
                      toType = 'SYMBOL', OrgDb = org.Hs.eg.db)


#--- fly id ---#
# auto_install('org.Dm.eg.db')
library(org.Dm.eg.db)
id = toTable(org.Dm.egSYMBOL) %>% dplyr::pull(1) %>% sample(20)
id = transId(id, trans_to = 'symbol',org='fly', return_dat = F)
transId(id, trans_to = 'ens',org='fly', return_dat = T)
