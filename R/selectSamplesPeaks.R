## select peaks for all samples
# INPUT:
# samples => samples as output of the function samplesSelection
# channel => channel to be considered, default="green"
# metric => metric to be used, selct among "median", "mean", "max" or "AUC", default is median
# Nchannel => channel to be used to normalise data, default is NA
# baseThr => threshold on the baseline used in order to define what is a peak, default is 0.01
# minLength => minimum length of a plug/droplet in number of data points, default is 10 (data points, not seconds)
# discartPeaks => select of to discart first and/or last peak ("first" discart the first, "last" discart the last, "both" discart both), default is NA
# discartPeaksPerc => select the percentage of peaks to discart if discartPeaks is defined. Default is 1.
# OUTPUT:
# samplesPeaks  => list with, in each position, a data.frame with max, mean, sd (standard deviation), start and length of all peaks

selectSamplesPeaks <- function(samples, channel="green", metric="median", Nchannel=NA, baseThr=0.01, minLength=10, discartPeaks=NA, discartPeaksPerc=1, Cchannel=NA){
  
  peaksALLsamples <- lapply(samples, peaksSelection, channel=channel, metric=metric, Nchannel=Nchannel, baseThr=baseThr, minLength=minLength, discartPeaks=discartPeaks, discartPeaksPerc= discartPeaksPerc, Cchannel=Cchannel)

	return(peaksALLsamples)    
}
