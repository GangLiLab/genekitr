rm(list=ls())
library(dplyr)
library(AnnoGenes)
data(geneList, package="DOSE")
ids = names(geneList)[1:100]
ensid = transId(ids, trans_to = 'ensembl',org = 'human')

gkeg <- genKEGG(ensid, org = 'hg',use_symbol = F)
head(gkeg)
plotEnrichDot(gkeg, xlab_type =  'FoldEnrich',
              show_item = 15, main_text_size = 14,legend_text_size = 10,
              low_color = 'red', high_color = 'blue',
              xleft = 0, font_type = 'Arial', remove_grid = T,
              wrap_width = 30,border_thick = 3 )

# if we use symbol gene as input:
smb = transId(ids,org = 'human',trans_to = 'symbol')
keg_raw <- genKEGG(smb, org = 'hs', use_symbol = F)
keg_readable <- genKEGG(smb, org = 'hs', use_symbol = T)


