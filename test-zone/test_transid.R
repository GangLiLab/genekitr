library(AnnoGenes)
data(geneList, package = 'DOSE')
ids = names(geneList)[1:10]
transId(id = ids, trans_to = 'symbol',org='hs', return_dat = T)
transId(id = ids, trans_to = 'ens',org='human', return_dat = F)

# 如果选择物种不对，会提示报错
transId(id = ids, trans_to = 'sym',org='human', return_dat = F)



