# data("allData", package="BraDiPluS")
load("/Users/eduati/BraDiPluS/data/allData_v2.RData")
setwd("/Users/eduati/BraDiPluS/fig/")

# compute and plot z-score for all cell line and patient samples 
res<-computeStatistics(allData, controlName="FS + FS", thPval=0.1, thSD=1.5, subsample=F, saveFiles=T)

# compute and plot z-score for patient #3
# exampleData <- list(patient3=allData$patient3)
BxPC3<-allData[2]
res<-computeStatistics(BxPC3, controlName="FS + FS", thPval=0.1, thSD=NA, subsample=F, saveFiles=T, showLabels=T)
