---
layout: page
title: Usage
---

## Preprocessing

This documentation is based on the example script (*/examples/workflow_preprocessing.R*) and data files (*data/BxPC3_data.RData*" and "*data/sampleNames.RData*") provided in the package.

For more detailed information look at the documentation for each funtion (typing "?nameOfFunction" in R).


For this example we will use screening data for the pancreatic cancer cell line BxPC3. Using the microfluidics platform cells were excapsulated in plugs together with drugs (alone or in combination) and a Caspase-3 fluorescent substrate, which becomes green fluorescent when cleaved by active Caspase-3. An automatially controlled Braille Display was used to generate a sequence of 62 experimental condition, each with 12 plugs (20 for the control condition without drugs), which can be considered as replicates. Experiemntal conditions (which will be called *samples* hereafter) were separated by blue fluorescent plugs used as barcode. An orange fluorescent dye was added to the cell suspension to verify the proper mixture of the components in each plug. Data are recorded using a Labview program developed by R. Utharala which can record three channels (*green*, *blue*, *orange*). During the measurement plugs are flushed at a constant flow rate and are therefore recorded as peaks in the corresponding channels. For more information of the experiment please refer to [the paper](http://biorxiv.org/content/early/2016/12/14/093906).


#### Load the data

You can load the data provided with the package as follow:

```R
data("BxPC3_data", package="BraDiPluS")
data("sampleNames", package="BraDiPluS")
```

The first line loads a dataframe called *MyData* with 4 columns corrisponging to measures in the *green*, *orange* and *blue* channels and acquisition *time* respectively. The second line loads a vector with corresponding sample names.

*MyData* includes one full cycle (will be called also *run*), which is a sequence of the 62 expeirmental conditions mentioned before. Normally the .txt file produced by the Labview program includes more than one sequence of the 62 samples, and these different runs can be separated using the function *cutData* and specifying start and end time of each run.

#### Visualise the data

To plot the data you can use the following function:

```R
plotData(data=MyData, channels=c("blue", "orange", "green"))
```

Please check the documentation ?plotData to see additional attributes which can be passed to the function to plot only a specific time range or to highlight info about peaks and samples which can be automatically detected with the functions described in the next sections.

#### Saparate samples and identify replicates for each sample

Different samples can be separated based on the peaks in the barcoding channel (which is the *blue* in our data):

```R
res <- samplesSelection(data=MyData, BCchannel="blue",BCthr=0.01,
	BCminLength=100, distThr=16, plotMyData=T, barcodePos="before")
samples<-res$samples
names(samples)<-sampleNames
```

Briefly, the samplesSelection function works as follows: 1. detect the barcoding peaks (default blue channel) and that is done considering the blue signal that goes above a defined threshold (set by argument *BCthr*) and with a minimum width (defined by *BCminLength* - which unit is number of data points). Then it looks at the space between selected blue peaks and consider as samples those with a distance higher than a threshold (defined by *distThr* - which unit is time). Arguments shoud be adjusted to properly select the samples. Selected samples can be also plotted using the *plotData* function.

```R
plotData(data=MyData, channels=c("blue", "orange", "green"), samples = samples)
```

Samples can then be processed to select the sequence of peaks in the *green* channel representig the replicates. In this step we can also discart the first peak, to remove potential cross contamination with previous sample, and short peaks, which correspond to unfused plugs.

```R
samplesPeaks <- selectSamplesPeaks(samples, channel="green", metric="median",
	baseThr=0.01, minLength=350, discartPeaks="first", discartPeaksPerc=5)
```

#### Quality assessment

After these steps, we considered the distribution of the intensity of the orange peaks across all samples and discarded the extreme values (i.e. the outliers). Where Q1 is the 25th percentile, Q3 is the 75th percentile and IQR is the interquantile range (Q3-Q1), outliers were defined as values lower than Q1 - 1.5*IQR, or higher than Q3+1.5*IQR. Before running the *qualityAssessment* function, samples with particularly high or low overall orange signal (e.g. due to cells clogging or malfunctioning inlet) were discarded:

```R
samplesPeaks[[indexSampleToRemove]]<-setNames(data.frame(matrix(ncol = ncol(samplesPeaks[[indexSampleToRemove]]), nrow = 0)), colnames(samplesPeaks[[indexSampleToRemove]]))
```

The function can than be run as:

```R
runs<-list(run1=samplesPeaks)
runs.qa<-qualityAssessment(runs=runs)
```

If more that one run (i.e. a cycle with all different experimental conditions) is available, and the corresponding peaks are called, for example *samplesPeaks1*, *samplesPeaks2*, *samplesPeaks3*, the function can be called as follows:

```R
runs<-list(run1=samplesPeaks1, run2=samplesPeaks2, run3=samplesPeaks3)
runs.qa<-qualityAssessment(runs=runs)
```

## Analysis

This part can be run together for all patient/cell line or separately. We will use as example the data provided with the package, that can be load as:

```R
data("allData", package="BraDiPluS")
```

Data are structured as a list 6 elements, one for each patient/cell line (type *names(allData)* to check the correspondig names). For example data for cell line BxPC3 can be accessed typing:

```R
BxPC3data <- allData$BxPC3
```

For each cell line/patient data are also structured as a list with one run (i.e. full repetition of all 62 experimetal conditions) for each element of the list. Each run contains the peaks selected using functions *selectSamplesPeaks* and *qualityAssessment*.

For example, we have 2 runs for patient #3 which can be accessed as:

```R
p3data_run1<-allData$patient3$run_1
p3data_run2<-allData$patient3$run_2
```

All runs for all patient/cell lines must have the same number of elements (which should correspond to the total number of experimental conditions, 62 in our case).


We will consider patient #3 to illustrate the usage of function *computeStatistics*.

```R
exampleData <- list(patient2=allData$patient2)
setwd("path/to/fuguresFolder")
res<-computeStatistics(exampleData, controlName="FS + FS", thPval=0.1, thSD=1.5, subsample=F, saveFiles=T)
```

This will produce two figures called "boxplot.pdf" and "volcanoplot.pdf" 


![boxplot patient 3](https://github.com/saezlab/BraDiPluS/blob/gh-pages/public/fig/boxplot_p3.jpg?raw=true)

