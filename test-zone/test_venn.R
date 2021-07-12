library(AnnoGenes)
library(dplyr)
set1 <- paste(rep("gene" , 100) , sample(c(1:1000) , 100 , replace=F) , sep="")
set2 <- paste(rep("gene" , 100) , sample(c(1:1000) , 100 , replace=F) , sep="")
set3 <- paste(rep("gene" , 100) , sample(c(1:1000) , 100 , replace=F) , sep="")
set4 <- paste(rep("gene" , 100) , sample(c(1:1000) , 100 , replace=F) , sep="")
set5 <- paste(rep("gene" , 100) , sample(c(1:1000) , 100 , replace=F) , sep="")

two_gene_list = list(gset1 = set1, gset2 = set2)
sm_gene_list = list(gset1 = set1, gset2 = set2, gset3 = set3)
la_gene_list = list(gset1 = set1, gset2 = set2, gset3 = set3, gset4 = set4, gset5 = set5 )

p1 = plotVenn(two_gene_list, alpha_degree = 1, border_thick = 0)

p2= plotVenn(sm_gene_list,alpha_degree = .3, border_thick = 1)
p3 = plotVenn(sm_gene_list,text_size = 2,alpha_degree = 1,remove_grid = T, color = ggsci::pal_lancet()(3))
p4 = plotVenn(la_gene_list,text_size = 10, border_thick = 2,remove_grid = T)

library(patchwork)
(p1+p2+p3)/p4
