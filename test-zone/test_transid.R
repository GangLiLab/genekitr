rm(list = ls())
library(AnnoGenes)
#--- huamn id ---#
data(geneList, package = 'DOSE')
id = names(geneList)[1:5]
id
transId(id, trans_to = 'symbol',org='hs')
transId(id, trans_to = 'uni',org='human')
transId(id, trans_to = 'ens',org='hg')

fake_id = c(id,'23215326','1','2','344263475','45')
transId(fake_id, trans_to = 'sym',org='human')

# 如果选择物种不对，会提示报错
transId(id, trans_to = 'sym',org='mouse')

# transId vs bitr
fake_id
res1 = transId(fake_id, trans_to = 'sym',org='human')
res2 = clusterProfiler::bitr(fake_id, fromType = 'ENTREZID',
                      toType = 'SYMBOL', OrgDb = org.Hs.eg.db)
class(res1); class(res2)
length(res1); nrow(res2)
res1
res2

#--- fly id ---#
# auto_install('org.Dm.eg.db')
library(org.Dm.eg.db)
id = toTable(org.Dm.egSYMBOL) %>% dplyr::pull(1) %>% sample(20)
transId(id, trans_to = 'ens',org='fly')
transId(id, trans_to = 'symbol',org='fly')



