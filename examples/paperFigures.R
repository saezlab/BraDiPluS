library(BraDiPluS)
# data("allData", package="BraDiPluS")
load("/Users/eduati/BraDiPluS/data/allData_v2.RData")
setwd("/Users/eduati/BraDiPluS/fig/")

###
## only for cell lines
cellLines<-allData[c(1,2)]

# compute and plot z-score for all cell line and patient samples 
res_stats<-computeStatistics(allData = cellLines, controlName="FS + FS", thPval=0.1, thSD=NA, subsample=F, saveFiles=F)

# compare conditions between samples 
# res<-compareSamples_volcanoplot(allData = cellLines, controlName="FS + FS", thPval=0.1, thEffectSize=0.3)
res<-compareSamples_barplot(allData = cellLines, controlName="FS + FS", thPval=0.05, res_stats=res_stats)



