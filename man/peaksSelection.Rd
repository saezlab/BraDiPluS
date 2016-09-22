\name{peaksSelection}
\alias{peaksSelection}

\title{
  Selection of fluorescence peaks
}

\description{
  This function allows to select peaks for the selected channel.
}

\usage{
peaksSelection(data, channel="green", metric="median", Nchannel=NA, baseThr=0.01, minLength=10, discartPeaks=NA, discartPeaksPerc=1, Cchannel=NA)
}


\arguments{
  \item{data}{data to be analyzed (formatted as data.frame with 4 columns: green, orange, blue, time)}
  \item{channel}{channel to be considered, default="green"}
  \item{metric}{metric to be used, selct among "median", "mean", "max" or "AUC", default is median}
  \item{Nchannel}{channel to be used to normalise data, default is NA}
  \item{baseThr}{threshold on the baseline used in order to define what is a peak, default is 0.01}
  \item{minLength}{minimum length of a plug/droplet in number of data points, default is 10 (data points, not seconds)}
  \item{discartPeaks}{select of to discart first and/or last peak ("first" discart the first, "last" discart the last, "both" discart both), default is NA}
  \item{discartPeaksPerc}{select the percentage of peaks to discart if discartPeaks is defined. Default is 1}
  \item{Cchannel}{if one channel was used for control (e.g. orange dye with cells), it can be used to discart peaks for which there was not a proper mixture, i.e. outliers in the control channel}
}
\details{
This function reads the data produced with labview in the 3 channels (green, orange and blue) and select the peaks corresponding to the droplets in the selected channel
}

\value{
This function returns a data.frame with one peak for each row and 9 columns:
\item{green}{value of the peak in the green channel}
\item{orange}{value of the peak in the orange channel}
\item{blue}{value of the peak in the blue channel}
\item{norm}{value of the selected channel normalized by the value of the Nchannel, if the Nchannel is provided (otherwise the value is set to 0)}
\item{start}{starting point of the peak}
\item{end}{final point of the peak}
\item{length}{length of the peak}
}

\author{
    F.Eduati
}

\seealso{
\link{plotData}
}

\examples{

data(BxPC3_data,package="BraDiPluS")
peaks <- peaksSelection(data=MyData, channel="blue")

}

