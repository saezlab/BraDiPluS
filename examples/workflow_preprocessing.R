data("BxPC3_data", package="BraDiPluS")
data("sampleNames", package="BraDiPluS")

# plot the data
plotData(data=MyData, channels=c("blue", "orange", "green"))

# separate the different samples (separated by barcode) 
res <- samplesSelection(data=MyData, BCchannel="blue",BCthr=0.01, BCminLength=100, distThr=16, plotMyData=F, barcodePos="before")
samples<-res$samples
names(samples)<-sampleNames

plotData(data=MyData, channels=c("blue", "orange", "green"), samples = samples)

# select the peaks for each sample
samplesPeaks <- selectSamplesPeaks(samples, channel="green", metric="median", baseThr=0.01, minLength=350, discartPeaks="first", discartPeaksPerc=5)

# remove outliers based on orange channel
runs<-list(run1=samplesPeaks)
runs.qa<-qualityAssessment(runs=runs)

# look at the median sample values
allData<-do.call(cbind, lapply(runs.qa, function(myRun){
  sapply(sapply(myRun, get, x="green"), median)
}))

# compute z-score
allData_scale<-apply(allData,2,scale)
