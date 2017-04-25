data("allData_v2", package="BraDiPluS")
setwd("/Users/eduati/BraDiPluS/fig/")

# compute and plot z-score for all cell line and patient samples 
res<-computeStatistics(allData, controlName="FS + FS", thPval=0.1, thSD=1.5, subsample=F, saveFiles=T)

# compute and plot z-score for patient #3
exampleData <- list(patient3=allData$patient3)
res<-computeStatistics(exampleData, controlName="FS + FS", thPval=0.1, thSD=NA, subsample=F, saveFiles=T, showLabels=T)
