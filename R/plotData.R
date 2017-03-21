#
#  This file is part of the `BraDiPluS` R package
#
#  Copyright (c) 2016 EMBL-EBI
#
#  File author(s): Federica Eduati (federica.eduati@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://github.com/saezlab/BraDiPluS
# --------------------------------------------------------
#
#' Plot the data and the selected peaks.
#' 
#' \code{plotData} allows to visualize the data.
#' 
#' This function is used to plot the desired range of data of one or more channels. It can be
#' used also to plot the peaks selected on one channel using the function \code{\link{peaksSelection}}.
#' 
#' @param data input data (formatted as data.frame with 4 columns: green, orange, blue, time),
#' imported from the .txt file saved by the in-house LabVIEW progam
#' @param channels a vector with the name of the channels we want to plot (e.g. channels=c("green", "blue")
#' to plot the green and the blue channel)
#' @param startTime initial time point to be plot. Default is NA, plots from the first time point in the data
#' @param ymin lower y limit. Default is NA, plots from the first time point in the data
#' @param ymax upper y limit. Default is NA, plots from the first time point in the data
#' @param peaks peaks computed using the function \code{\link{peaksSelection}}
#' @param samples samples selected using the function \code{\link{samplesSelection}}
#' @param ytext position in the y axes of the sample names/numbers if samples argument is specified
#' @param savePlot if FALSE, the plot is not saved; if TRUE plot is saved with default file name; if string,
#' plot is saved with given filename. Plot is saved in .png format
#' @seealso \code{\link{peaksSelection}} to select peaks in a specific channel, \code{\link{samplesSelection}} to separate
#' data in samples
#' @examples 
#' data(BxPC3_data, package="BraDiPluS")
#' plotData(data=MyData, channels=c("blue", "green", "orange"))
#' 
#' res <- samplesSelection(data=MyData, BCchannel="blue",BCthr=0.01, BCminLength=100, distThr=16, plotMyData=F, barcodePos="before")
#' samples<-res$samples
#' plotData(data=MyData, channels=c("blue", "orange", "green"), samples = samples)
#' 
#' BCpeaks <- peaksSelection(MyData, channel="blue")
#' plotData(data=MyData, channels=c("blue", "orange", "green"), peaks = BCpeaks)
#' @export

plotData <- function(data, channels, startTime=NA, endTime=NA, ymin=NA, ymax=NA, peaks=NA, samples=NA, ytext=(-0.005), savePlot=FALSE){
	

  if (savePlot!=FALSE){
  	if (savePlot==TRUE){
  		filename<-paste("plotData", deparse(substitute(data)), sep="_")
  		filename<-paste(filename, "png", sep=".")
  	}else{
  		filename<-savePlot
  	}
  	ppi <- 300; png(filename, width=16*ppi, height=9*ppi, res=ppi)
  }
	
  if (nrow(data)==0){
    plot(1, type="n", axes=F, xlab="", ylab="")
    text(1,1,"NO DATA", cex=3)
  }else{
    if (is.na(startTime)){startTime<-data$time[1]}
    if (is.na(endTime)){endTime<-data$time[nrow(data)]}
    
    ix1 <- which(data$time>=startTime)
    ix2 <- which(data$time>=endTime)
    
    col <- rep(NA, length(channels))
    names(col) <- channels
    if ("green" %in% channels) {col["green"] <- "#31a354"}
    if ("orange" %in% channels) {col["orange"] <- "#e6550d"}
    if ("blue" %in% channels) {col["blue"] <- "#3182bd"}  
    
    data <- data[ix1[1]:ix2[1],]
    if (is.na(ymax)){ymax <- max(data[,channels]) + abs(max(data[,channels])/100)} 
	if (is.na(ymin)){ymin <- min(data[,channels]) - abs(min(data[,channels])/100)} 
    
    plot(data$time, as.matrix(data[channels[1]]), ylim=c(ymin, ymax), col=col[1], type="l", xlab="time (sec)", ylab="fluorescence", main="", cex.lab=1.5)
    if (length(channels) > 1){
      for (i in 2:length(channels))
        lines(data$time, as.matrix(data[channels[i]]), col=col[i], type="l", xlab="time", ylab="fluorescence", main="")
    }
    
    if (all(!is.na(peaks))){
      if (nrow(peaks)>0){
        tmp_start <- which(peaks$start>startTime)[1]
        tmp_end <- which(peaks$start>endTime)[1]
        if (!is.na(tmp_end)){
          tmp_end <- tmp_end -1
        } else{
          tmp_end <- nrow(peaks)
        }
        peaks <- peaks[tmp_start:tmp_end,]
        
        pch <- rep(NA, length(channels))
        names(pch) <- channels
        if ("green" %in% channels) {pch["green"] <- 21}
        if ("orange" %in% channels) {pch["orange"] <- 22}
        if ("blue" %in% channels) {pch["blue"] <- 23}  
        points(peaks$start, as.matrix(peaks[channels[1]]), bg=col[1], pch=pch[1])
        if (length(channels) > 1){
          for (i in 2:length(channels))
            points(peaks$start, as.matrix(peaks[channels[i]]), bg=col[i], pch=pch[i])
        }
      }
      
      
      #points(peaks$start, peaks$median, col="#de2d26")
      #text(peaks$start, peaks$median+200, col="#de2d26")
    }
    
    if (all(!is.na(samples))){
      ytext2<-(ytext+1.5*(ytext))
      plotSamples<-function(x){
        if (length(x)>0){
          lines(x=c(x[1], x[length(x)]), y=c(ytext, ytext))
        }
        
      }
      tmp<-lapply(lapply(samples, get, x="time"), plotSamples)
      #     text(lapply(lapply(samples, get, x="time"), function(x) (x[1]+(x[length(x)]-x[1])/2)), y=ytext2, labels=seq(1, length(samples)), cex=0.5)
      text(lapply(lapply(samples, get, x="time"), function(x) (x[1]+(x[length(x)]-x[1])/2)), y=ytext2, labels=names(samples), cex=0.5)    
    }
    
  }
  
  
  if (savePlot!=FALSE){dev.off()}

}