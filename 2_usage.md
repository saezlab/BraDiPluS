---
layout: page
title: Usage
---

This documentation is based on the example script (/examples/workflow_preprocessing.R) and data files ("data/BxPC3_data.RData" and "data/sampleNames.RData") provided in the package.


Load the example data provided in the package

```R
data("BxPC3_data", package="BraDiPluS")
data("sampleNames", package="BraDiPluS")
```

plot the data

```R
plotData(data=MyData, channels=c("blue", "orange", "green"))
```

separate the different samples (separated by barcode) 

```R
res <- samplesSelection(data=MyData, BCchannel="blue",BCthr=0.01, BCminLength=100, distThr=16, plotMyData=F, barcodePos="before")
samples<-res$samples
names(samples)<-sampleNames

plotData(data=MyData, channels=c("blue", "orange", "green"), samples = samples)
```

select the peaks for each sample

```R
samplesPeaks <- selectSamplesPeaks(samples, channel="green", metric="median",
							baseThr=0.01, minLength=350, discartPeaks="first", discartPeaksPerc=5)
```

remove outliers based on orange channel

```R
runs<-list(run1=samplesPeaks)
runs.qa<-qualityAssessment(runs=runs)
```

look at the median sample values

```R
allData<-do.call(cbind, lapply(runs.qa, function(myRun){
sapply(sapply(myRun, get, x="green"), median)
}))
```

compute z-score

```R
allData_scale<-apply(allData,2,scale)
```
