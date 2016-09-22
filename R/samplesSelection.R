## separate samples
# INPUT:
# data => all data to be analyzed (formatted as data.frame with at least 2 columns: time and blue (blue channel))
# BCchannel =>channel used for barcoding, default is "blue"
# BCthr <- threshold used in order to define what is a peak. Default is 0.01
# distThr <- minimum distance between peaks to be considered a sample
# barcodePos <- specify if the barcode is before or after the sample ("before" or "after"), default is "after"
# OUTPUT:
# samples => list containing a sample (data.frame with arguments frame and FL) in each position
# barcode => corresponding barcode for each sample 


samplesSelection <- function(data, channel="green", BCchannel="blue", BCthr=0.01, distThr=100, plotMyData=FALSE, BCminLength=10, barcodePos="after"){
  
  # slect peaks in the barcoding channel
  peaks <- peaksSelection(data, channel=BCchannel, baseThr=BCthr, minLength=BCminLength)
  
  if (plotMyData==TRUE){plotData(data, channels=c(BCchannel, channel), peaks=peaks)}
  
  # distance of each peak to the previous one
  peaksDist<-c(0, peaks$start[2:length(peaks$start)] - peaks$end[1:length(peaks$start)-1])
  # these are the indices of the first peak of each barcode (except for the first BC) 
  ix <- c(which(peaksDist>distThr))
  
  if (length(ix)==0) stop(paste("Error: check the argument 'distThr', your maximum distance between peaks is", max(peaksDist)))
  
  
  ####
  # define beginning and end of the BC samples 
  BCstart.time<-c(peaks$start[1], peaks$start[ix])
  #abline(v=BCstart.time, col="red")
  BCstart.ix<-sapply(BCstart.time,function (x) which(data$time>=x)[1]) - 10
  BCstart.ix[BCstart.ix<=1]<-1
    
  BCend.time<-c(peaks$end[ix-1], peaks$end[nrow(peaks)])
  #abline(v=BCend.time, col="blue")
  BCend.ix<-sapply(BCend.time,function (x) which(data$time>=x)[1]) + 10
  BCend.ix[BCstart.ix>=nrow(data)]<-nrow(data)
  
  ####
  # define beginning and end of the real sample
  RSstart.time<-c(peaks$end[ix-1])
  #abline(v=RSstart.time, col="blue")
  RSstart.ix<-sapply(RSstart.time,function (x) which(data$time>=x)[1]) - 10
  RSstart.ix[RSstart.ix<=1]<-1
  
  RSend.time<-c(peaks$start[ix])
  #abline(v=RSend.time, col="red")
  RSend.ix<-sapply(RSend.time,function (x) which(data$time>=x)[1]) + 10
  RSend.ix[RSend.ix>=nrow(data)]<-nrow(data)

  
  samples <- list()
  for (i in 1:length(RSstart.time)){
    samples[[i]]<-data[RSstart.ix[i]:RSend.ix[i],]
    names(samples)[i] <- paste("sample", i, sep="_")
    
    if (plotMyData==TRUE){
      # colour only border
      tmp <- data[,c(channel, BCchannel)]
      rect(RSstart.time[i],0,RSend.time[i],max(tmp), border="#31a354") 
      text(x=(RSstart.time[i]+(RSend.time[i]-RSstart.time[i])/2), y=max(tmp), labels=i, cex=0.8)
    }
  }
  
  barcode <- list()
  for (i in 1:length(BCstart.time)){
    barcode[[i]]<-data[BCstart.ix[i]:BCend.ix[i],]
    names(barcode)[i] <- paste("sample", i, sep="_")
    
#     if (plotMyData==TRUE){
#       # colour only border
#       tmp <- data[,c(BCchannel)]
#       rect(BCstart.time[i],0,BCend.time[i],max(tmp), border="#2F4F4F") 
#       text(x=(BCstart.time[i]+(BCend.time[i]-BCstart.time[i])/2), y=max(tmp), labels=i, cex=0.8)
#     }
    
  }
  
  
  # the barcoding starts at the end of the sample
  if (barcodePos=="after"){
    barcode<-barcode[2:length(barcode)]
  }else if (barcodePos=="before"){
    barcode<-barcode[1:(length(barcode)-1)]
  }else{
    stop(cat("Please specify a valid argument for barcode, can be: \'before\' or \'after\'"))
  }
  
  
  
  return(list(samples=samples, barcode=barcode))
}
  