BraDiPluS
=========
BraDiPluS - Braille Display Plugs Data Analysis - is an R package developed to analyse the data of plugs produced using in-house combinatorial microfluidics platform and adquired using an in-house LabVIEW program allowing the detection of three channels:
* _blue_: sample barcodes
* _green_: Caspase-3 activity assay
* _orange_: marker dye to monitor mixing of reagents 


## How to install it

Open R and, from the folder with the package, type:

```
install.packages('BraDiPluS',repos=NULL,type='src')
```

## Workflow example

Load the example data provided in the package

```
data("BxPC3_data", package="BraDiPluS")
data("sampleNames", package="BraDiPluS")
```

plot the data

```
plotData(data=MyData, channels=c("blue", "orange", "green"))
```

separate the different samples (separated by barcode) 

```
res <- samplesSelection(data=MyData, BCchannel="blue",BCthr=0.01, BCminLength=100, distThr=16, plotMyData=F, barcodePos="before")
samples<-res$samples
names(samples)<-sampleNames

plotData(data=MyData, channels=c("blue", "orange", "green"), samples = samples)
```

select the peaks for each sample
```
samplesPeaks <- selectSamplesPeaks(samples, channel="green", metric="median", baseThr=0.01, minLength=350, discartPeaks="first", discartPeaksPerc=5)
```

remove outliers based on orange channel
```
runs<-list(run1=samplesPeaks)
runs.qa<-qualityAssessment(runs=runs)
```

look at the median sample values
```
allData<-do.call(cbind, lapply(runs.qa, function(myRun){
sapply(sapply(myRun, get, x="green"), median)
}))
```

compute z-score
```
allData_scale<-apply(allData,2,scale)
```


![alt text](https://s-media-cache-ak0.pinimg.com/564x/50/0a/43/500a43906dea54dbff4d81d9ab7dfb23.jpg "BraDiPluS" width="10")