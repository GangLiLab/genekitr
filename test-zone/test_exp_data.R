devtools::install_github("xjsun1221/tinyarray",upgrade = F)
library(tinyarray)
?get_deg()

BiocManager::install('airway')
library(airway)
data(airway)
exprSet=assay(airway)
group_list = airway$dex
library(DESeq2)
res = des(exprSet,group_list)
deg = res[[1]]
