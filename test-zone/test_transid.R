rm(list = ls())
library(genekitr)
#--- huamn id ---#
data(geneList, package = 'genekitr')
id = names(geneList)[1:100]
id
transId(id, trans_to = 'symbol',org='hs')
transId(id, trans_to = 'uni',org='human')
transId(id, trans_to = 'ens',org='hg')

fake_id = c(id,'23215326','1','2','344263475','45')
transId(fake_id, trans_to = 'sym',org='human')

# 如果选择物种不对，会提示报错
transId(id, trans_to = 'sym',org='mouse')

# transId vs bitr
fake_id = c("DEL11P13" , "TRNAV-CAC", "MMD2" ,     "HBD"  ,     "RNR1",
            "RNR2" ,     "TEC"  ,     "MEMO1" ,    "TP53"  ,    "BCC7")

res1 = transId(fake_id, trans_to = 'ens',org='human',unique = T)
res2 = clusterProfiler::bitr(fake_id, fromType = 'SYMBOL',
                      toType = 'ENSEMBL', OrgDb = org.Hs.eg.db)

res3 = transId(fake_id, trans_to = 'entrez',org='human',unique = F)
res4 = clusterProfiler::bitr(fake_id, fromType = 'SYMBOL',
                             toType = 'ENTREZID', OrgDb = org.Hs.eg.db)


#--- fly id ---#
# auto_install('org.Dm.eg.db')
library(org.Dm.eg.db)
id = toTable(org.Dm.egSYMBOL) %>% dplyr::pull(1) %>% sample(20)
transId(id, trans_to = 'ens',org='fly')
transId(id, trans_to = 'symbol',org='fly')
transId(id, trans_to = 'uni',org='fly')

#--- MasSpec uniprot to symbol ---#
load('./test-zone/phosphosite_test.Rdata')
uni1 = x$`Protein accession`
sym1 = x$`Gene name`

system.time({
  sym2 = transId(uni1,'sym','mm')
})
which(!sym2 %in% sym1)

system.time({
  sym3 = clusterProfiler::bitr(uni1,fromType = 'UNIPROT',toType = 'SYMBOL',OrgDb = org.Mm.eg.db)
})
length(sym2);dim(sym3) #get much more symbols for identical uniprot id

# exp: Q76KJ5 https://www.uniprot.org/uniprot/Q76KJ5
x[1524,]
sym2[1524]

# exp2: Q6P2L7  https://www.uniprot.org/uniprot/Q6P2L7
x[1355,]
sym2[1355]



