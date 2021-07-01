data(geneList, package="DOSE")
head(geneList)

source('./R/genGSEA.R')
source('./R/genInfo.R')
head(geneList)
genInfo(names(head(geneList)), org = 'hs') %>%
  dplyr::pull(symbols)

genGSEA(genelist = geneList, species = 'human',
        category = 'C1')
