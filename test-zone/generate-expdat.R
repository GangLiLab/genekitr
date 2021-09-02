library("airway")
data(airway)
airway
as.data.frame(colData(airway))
summary(colSums(assay(airway))/1e6)
metadata(rowRanges(airway))
