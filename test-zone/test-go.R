rm(list=ls())
library(dplyr)
library(ggplot2)
library(AnnoGenes)
data(geneList, package="DOSE")
id = names(geneList)[1:500]

transId(id, trans_to = 'symbol',org = 'human')

ego <- genGO(id, org = 'human',ont = 'mf',pvalueCutoff = 0.05,qvalueCutoff = 0.1 ,use_symbol = F)
head(ego)

test = rio::import('~/Downloads/GO_analysis_EGR1_05.xlsx')
p1 = plotEnrichDot(test, xlab_type =  'FoldEnrich', legend_by = 'qvalue',
              show_item = 15, main_text_size = 14,legend_text_size = 10,
              low_color = 'red', high_color = 'blue',
              xleft = 0, font_type = 'Arial', remove_grid = T,
              wrap_width = 30,border_thick = 3 )

p2=plotEnrichDot(ego, xlab_type =  'GeneRatio', legend_by = 'p.adjust',
                 show_item = 10, main_text_size = 14,legend_text_size = 10,
                 low_color = 'orange', high_color = 'green',
                 xleft = 0, font_type = 'Times New Roman', remove_grid = F,
                 wrap_width = NULL ,border_thick = 1)
library(patchwork)
p1+ p2






