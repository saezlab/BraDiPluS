\name{plotData}
\alias{plotData}
\title{
Plot data and selected peaks
}
\description{
This function is used to plot the desired range of data of one or more channels. It can be used also to plot the peaks selected on one channel using the function \link{peaksSelection}.
}
\usage{
plotData(data, channels, startTime=NA, endTime=NA, ymin=NA, ymax=NA, peaks=NA, samples=NA, ytext=(-0.005), savePlot=FALSE)
}
\arguments{
  \item{data}{
input data, must have at least a column time and one of the channels (green, orange and/or blue)
}
  \item{channels}{
a vector with the name of the channels we want to plot (e.g. channels=c("green", "blue") to plot the green and the blue channel)
}
  \item{startTime}{
initial time point to be plot. Default is NA, plots from the first time point in the data
}
  \item{ymin}{
lower y limit. Default is NA, plots from the first time point in the data
}
  \item{ymax}{
upper y limit. Default is NA, plots from the first time point in the data
}
  \item{peaks}{
peaks computed using the function peaksSelection
}
  \item{samples}{
samples selected using the function samplesSelection
}
  \item{ytext}{
 position in the y axes of the sample names/numbers if samples argument is specified 
}
   \item{savePlot}{
if FALSE, the plot is not saved; if TRUE plot is saved with default file name; if string, plot is saved with given filename. Plot is saved in .png format
}
}

\value{

}
\author{
F.Eduati
}
\seealso{
\link{peaksSelection}, \link{samplesSelection}
}
\examples{
data(BxPC3_data, package="BraDiPluS")
plotData(data=MyData, channels=c("blue", "green", "orange"))

res <- samplesSelection(data=MyData, BCchannel="blue",BCthr=0.01, BCminLength=100, distThr=16, plotMyData=F, barcodePos="before")
samples<-res$samples
plotData(data=MyData, channels=c("blue", "orange", "green"), samples = samples)

}
