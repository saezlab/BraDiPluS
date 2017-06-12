library(BraDiPluS)
# data("allData", package="BraDiPluS")
load("/Users/eduati/BraDiPluS/data/allData_v2.RData")
setwd("/Users/eduati/BraDiPluS/fig/")

res_stats<-computeStatistics(allData = allData, controlName="FS + FS", thPval=0.05, thSD=NA, subsample=T, saveFiles=F)


###
## only for cell lines
cellLines<-allData[c(1,2)]

# compute and plot z-score for all cell line and patient samples 
res_stats<-computeStatistics(allData = cellLines, controlName="FS + FS", thPval=0.05, thSD=NA, subsample=F, saveFiles=F)

AsPC1_OK<-subset(res_stats, ((patientCellLine=="cell line: AsPC1") & (threshold=="Sign")))
BxPC3_OK<-subset(res_stats, ((patientCellLine=="cell line: BxPC3") & (threshold=="Sign")))

nrow(AsPC1_OK)
nrow(BxPC3_OK)


# compare conditions between samples 
# res<-compareSamples_volcanoplot(allData = cellLines, controlName="FS + FS", thPval=0.1, thEffectSize=0.3)
res<-compareSamples_barplot(allData = cellLines, controlName="FS + FS", thPval=0.05, res_stats=res_stats)


###
## only for patients
patients<-allData[c(3,4,5,6)]

# compute and plot z-score for all cell line and patient samples 
res_stats<-computeStatistics(allData = patients, controlName="FS + FS", thPval=0.1, thSD=NA, subsample=F, saveFiles=F)

# compare conditions between samples 
# res<-compareSamples_volcanoplot(allData = cellLines, controlName="FS + FS", thPval=0.1, thEffectSize=0.3)
res<-compareSamples_barplot(allData = patients, controlName="FS + FS", thPval=0.05, res_stats=res_stats)

