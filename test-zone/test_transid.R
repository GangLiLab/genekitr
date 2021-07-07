library(AnnoGenes)
data(geneList, package = 'DOSE')
ids = names(geneList)[1:10]
transId(id = ids, trans_to = 'symbol',org='hs', return_dat = T)
transId(id = ids, trans_to = 'symbol',org='human', return_dat = F)
transId(id = ids, trans_to = 'ens',org='human', return_dat = F)

# 如果选择物种不对，会提示报错
transId(id = ids, trans_to = 'symbol',org='mouse', return_dat = F)

ids =c('ENSG00000260179','ENSG00000234396','ENSG00000225972')
human_gtf %>% filter(ensembl %in% ids) %>%
  select(ensembl,symbol)

x = transId(id = human_gtf$ensembl, trans_to = 'entrez',org='hs', return_dat = T)
