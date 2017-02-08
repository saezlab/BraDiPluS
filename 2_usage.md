---
layout: page
title: Usage
---

This documentation is based on the example script (/examples/workflow_preprocessing.R) and data files ("data/BxPC3_data.RData" and "data/sampleNames.RData") provided in the package.

For more detailed information look at the documentation for each funtion (typing "?nameOfFunction" in R).


For this example we will use screening data for the pancreatic cancer cell line BxPC3. Using the microfluidics platform cells were excapsulated in plugs together with drugs (alone or in combination) and a Caspase-3 fluorescent substrate, which becomes green fluorescent when cleaved by active Caspase-3. An automatially controlled Braille Display was used to generate a sequence of 62 experimental condition, each with 12 plugs (20 for the control condition without drugs), which can be considered as replicates. Experiemntal conditions (which will be called *samples* hereafter) were separated by blue fluorescent plugs used as barcode. An orange fluorescent dye was added to the cell suspension to verify the proper mixture of the components in each plug. Data are recorded using a Labview program developed by R. Utharala which can record three channels (*green*, *blue*, *orange*). During the measurement plugs are flushed at a constant flow rate and are therefore recorded as peaks in the corresponding channels. For more information of the experiment please refer to [the paper](http://biorxiv.org/content/early/2016/12/14/093906).


#### Load the data

You can load the data provided with the package as follow:

```R
data("BxPC3_data", package="BraDiPluS")
data("sampleNames", package="BraDiPluS")
```

The first line loads a dataframe called *MyData* with 4 columns corrisponging to measures in the *green*, *orange* and *blue* channels and acquisition *time* respectively. The second line loads a vector with corresponding sample names.

#### Visualise the data

To plot the data you can use the following function:

```R
plotData(data=MyData, channels=c("blue", "orange", "green"))
```

![plot of all data](https://github.com/saezlab/BraDiPluS/blob/gh-pages/public/fig/allData.png)

Please check the documentation ?plotData to see additional attributes which can be passed to the function to plot only a specific time range or to highlight info about peaks and samples which can be automatically detected with the functions described in the next sections.

#### Identify replicates for each sample

Different samples can be separated based on the peaks in the barcoding channel (which is the *blue* in our data):

```R
res <- samplesSelection(data=MyData, BCchannel="blue",BCthr=0.01,
	BCminLength=100, distThr=16, plotMyData=F, barcodePos="before")
samples<-res$samples
names(samples)<-sampleNames

plotData(data=MyData, channels=c("blue", "orange", "green"), samples = samples)
```

Samples can then be processed to select the sequence of peaks in the *green* channel representig the replicates. In this step we can also discart the first peak, to remove potential cross contamination with previous sample, and short peaks, which correspond to unfused plugs.

```R
samplesPeaks <- selectSamplesPeaks(samples, channel="green", metric="median",
	baseThr=0.01, minLength=350, discartPeaks="first", discartPeaksPerc=5)
```

#### Quality assessment

After these steps, we considered the distribution of the intensity of the orange peaks across all samples and discarded the extreme values (i.e. the outliers). Where Q1 is the 25th percentile, Q3 is the 75th percentile and IQR
600 is the interquantile range (Q3-Q1), outliers were defined as values lower than Q1 - 1.5*IQR, or higher than Q3+1.5*IQR. These strict rules were applied to guarantee higher quality of the data used for the analysis described in the paper.

```R
runs<-list(run1=samplesPeaks)
runs.qa<-qualityAssessment(runs=runs)
```

For each sample, median value and z-score can be computed as follow:

```R
allData<-do.call(cbind, lapply(runs.qa, function(myRun){
sapply(sapply(myRun, get, x="green"), median)
}))

allData_zscore<-apply(allData,2,scale)
```

