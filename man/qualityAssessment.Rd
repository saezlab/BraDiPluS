\name{qualityAssessment}
\alias{qualityAssessment}

\title{
  Remove outlier peaks.
}

\description{
  This function allows to remove the peaks that are outliers based on the orange channel, meaning that the mixture of the component was not correct.
}

\usage{
qualityAssessment(runs)
}


\arguments{
  \item{runs}{List with one element for each run. For each run the list element correspond to the peaks selected for the different samples using \link{selectSamplesPeaks}}
 }
\details{
This function select remove the outlier peaks for all the samples previously identified using \link{selectSamplesPeaks}.
}

\value{
This function returns a list with the same structure of the one used as input but without the outlier peaks.
}

\author{
    F.Eduati
}

\seealso{
\link{samplesSelection}, \link{samplesBoxPlot}, \link{samplesBarPlot}
}

\examples{

data(BxPC3_data,package="BraDiPluS")
res <- samplesSelection(data=MyData, BCchannel="blue",
						BCthr=0.01, distThr=300, plotMyData=TRUE)
samples <-res$samples

# select the peaks for each sample
samplesPeaks <- selectSamplesPeaks(samples, channel="green", metric="median", baseThr=0.01, minLength=350, discartPeaks="first", discartPeaksPerc=5)

# remove outliers based on orange channel
runs<-list(run1=samplesPeaks)
runs.qa<-qualityAssessment(runs=runs)
}

