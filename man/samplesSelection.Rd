\name{samplesSelection}
\alias{samplesSelection}

\title{
  Selection of of the samples using the barcoding channel to separate them
}

\description{
  This function allows to select the samples using the selected barcoding channel.
}

\usage{
samplesSelection(data, channel="green", BCchannel="blue",
				 BCthr=0.01, distThr=100, plotMyData=FALSE, BCminLength=10, barcodePos="after")
}


\arguments{
  \item{data}{data to be analyzed (formatted as data.frame with 4 columns: green, orange, blue, time)}
  \item{channel}{main channel to be shown, default="green"}
  \item{BCchannel}{channel used for barcoding, default is "blue"}
  \item{BCthr}{threshold used in order to define what is a peak in the barcoding channel. Default is 0.01}
  \item{distThr}{minimum duration of a sample (in the same unit used for time)}
  \item{plotMyData}{TRUE to show the plot of the data, FALSE not to show it}
  \item{BCminLength}{minimum length of a plug/droplet in number of data points for the barcoding channel, default is 10 (data points, not seconds}
  \item{barcodePos}{specify if the barcode is before or after the sample ("before" or "after"), default is "after"}
}
\details{
This function selects the samples and optionally shows them.
}

\value{
This function returns a list with two components:
\item{samples}{list containing a sample in each position. Each sample is a data.frame with the same structure of the input data}
\item{barcode}{corresponding barcode for each sample. Each barcode is a data.frame with the same structure of the input data}
}

\author{
    F.Eduati
}

\seealso{
\link{selectSamplesPeaks}
}

\examples{

data(BxPC3_data,package="BraDiPluS")
res <- samplesSelection(data=MyData, BCchannel="blue",BCthr=0.01
				 BCminLength=100, distThr=16, barcodePos="before")
samples <-res$samples
barcode <-res$barcode

}

