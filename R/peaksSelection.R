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
#' Select fluorescence peaks.
#' 
#' \code{peaksSelection} allows to identify peaks for the defined channel.
#' 
#' This function reads the data produced with labview in the 3 channels (green, orange and blue)
#' and selects the peaks corresponding to the plugs in the user defined channel. This is done by considering
#' as peaks the signal that goes above a threshold (defined by argument "baseThr") and stays above for a minimum number
#' of data points (defined by argument "minLength").
#' 
#' @param data data to be analyzed (formatted as data.frame with 4 columns: green, orange, blue, time)
#' @param channel channel to be considered, default="green"
#' @param metric metric to be used to define the value to assign to each peak, select among "median", "mean",
#' "max" or "AUC", default is "median"
#' @param Nchannel channel to be used to normalise data, default is NA (raccomended to not use it)
#' @param baseThr threshold on the baseline used in order to define what is a peak, default is 0.01
#' @param minLength minimum length of a plug/droplet in number of data points, default is 10 (note that unit is
#' number of data points, not seconds)
#' @param discartPeaks select of to discart first and/or last peak ("first" discart the first, "last" discart the last,
#' "both" discart both), default is NA
#' @param discartPeaksPerc select the percentage of peaks to discart if discartPeaks is defined. Default is 1
#' @return This function returns a data.frame with one peak for each row and 9 columns:
#' \describe{
#'   \item{green}{value of the peak in the green channel}
#'   \item{orange}{value of the peak in the orange channel}
#'   \item{blue}{value of the peak in the blue channel}
#'   \item{norm}{value of the selected channel normalized by the value of the Nchannel, if the Nchannel is provided (otherwise the value is set to 0)}
#'   \item{start}{starting point of the peak}
#'   \item{end}{final point of the peak}
#'   \item{length}{length of the peak}
#' }
#' @seealso \code{\link{plotData}}
#' @examples 
#' data(BxPC3_data,package="BraDiPluS")
#' peaks <- peaksSelection(data=MyData, channel="blue")
#' @export

peaksSelection <- function(data, channel="green", metric="median", Nchannel=NA, baseThr=0.01, minLength=10, discartPeaks=NA, discartPeaksPerc=1){
  
  # consider only data for the desired channel
  x<-as.matrix(data[,channel])
  
  if (length(x)==0){
    cat("CAREFUL: no peaks identified for the current sample - no data\n")
    theColNames<-c("green", "orange", "blue", "norm", "start", "end", "length")
    peaks <- setNames(data.frame(matrix(ncol = length(theColNames), nrow = 0)), theColNames)
  }else{
    # define when data are above a selected threshold
    thr <- (x>baseThr) #vector with TRUE for values above the baseline threshold and FALSE for values below
    # find the corresponding indices
    
    # to find the start points I check when the T is preceeded by a non T (i.e. a F)
    ix_T<-c(0, which(thr==TRUE))
    N<-length(ix_T)
    ix_T_diff<-ix_T[2:N]-ix_T[1:(N-1)]
    ix_start<-ix_T[(which(ix_T_diff!=1))+1]
    
    
    # to find the end points I check when the F is preceeded by a non F (i.e. a T)
    ix_F<-which(thr==FALSE)
    N<-length(ix_F)
    ix_F_diff<-ix_F[2:N]-ix_F[1:(N-1)]
    ix_end<-ix_F[(which(ix_F_diff!=1))+1]-1
    
    # if ix_start is longer than ix_end it means that the last peak starts but do not end so I remove it
    ix_start<-ix_start[1:length(ix_end)]
    
    # compute corresponding length
    pLength<-ix_end-ix_start+1
    
    # remove short peaks (below the defined threshold)
    ix_start <- ix_start[which(pLength > minLength)]
    ix_end <- ix_end[which(pLength > minLength)]
    pLength <- pLength[which(pLength > minLength)]
    
    
    # check if there are peaks identified 
    if (length(ix_start)==0){
      cat("CAREFUL: no peaks identified for the current sample\n")
      theColNames<-c("green", "orange", "blue", "norm", "start", "end", "length")
      peaks <- setNames(data.frame(matrix(ncol = length(theColNames), nrow = 0)), theColNames)
    }else{
      # select corresponding data and compute metric
      x_all<-data
      x_all$time<-NULL
      
      pAll<-NA*x_all[1:length(ix_start),]
      rownames(pAll)<-NULL
      
      if (metric!="AUC"){ 
        metric <- match.fun(metric)
        for (i in 1:length(ix_start)){
          tmp_all<- x_all[seq(ix_start[i], ix_end[i]), , drop = FALSE]
          pAll[i,]<-apply(tmp_all,2,metric) # store the computed metric for each channel
        }
        
      }else{ # AUC is computed as the sum of data points multiplied by the sampling time
        metric<-function(x, end.time, start.time){sum(x)*(end.time-start.time)}
        for (i in 1:length(ix_start)){
          tmp_all <- x_all[seq(ix_start[i], ix_end[i]), , drop = FALSE]
          pAll[i,]<-apply(tmp_all,2,metric, data$time[ix_end][i], data$time[ix_start][i])
        }
      }
      
      # store the corresponding time points
      pStart <- data$time[ix_start]
      pEnd <- data$time[ix_end]
      
      # if the normalization channel is provided, use it to normalize the data
      if (!is.na(Nchannel)){
        pNorm<-pAll[,channel]/pAll[,Nchannel]
      } else{
        pNorm<-rep(0, nrow(pAll))
      }
      
      # save relevant values 
      peaks <- data.frame(green=pAll$green, orange=pAll$orange, blue=pAll$blue, norm=pNorm, start=pStart, end=pEnd, length=pLength) 
      
      
      # and remove short peaks (below the defined threshold)
      #peaks <- peaks[which(peaks$length > minLength),]
      
      # discart first and/or last replicates because of contamination
      if (!is.na(discartPeaks)){
        ix1<-1
        ix2<-nrow(peaks)
        nRemove<-ceiling(ix2* discartPeaksPerc*0.01) # number of peaks to remove according to the defined discartPeaksPerc
        if (discartPeaks=="first"){
          if (ix2>=2){
            ix1<-1+nRemove
          }
          else{warning("there is only one peak selected, cannot discart first peak")}
        }
        if (discartPeaks=="last"){
          if (ix2>=2) {
            ix2<-(nrow(peaks)-nRemove)
          }
          else{warning("there is only one peak selected, cannot discart last peak")}
        }
        if (discartPeaks=="both"){
          if (ix2>=3) {
            ix1<-1+nRemove;
            ix2<-(nrow(peaks)-nRemove)
          }
          else{warning("there are less then 2 peaks selected, cannot discart first and last peak")}
        }
        
        if (ix1>ix2){ix2<-ix1;warning("if you are discarting peaks from both sides the percentage must be lower then 50")}
        
        peaks <- peaks[ix1:ix2,]
      }
      
      # green shows in orange, thus if green is the main channel and orange the mixture control channel
      # orange is corrected as follow orange_new <- orange - 0.45*green
      
      #   peaks$orangeOriginal<-peaks$orange
      #   if (!is.na(Cchannel)){
      #   	if (channel=="green" & Cchannel=="orange"){
      #   	peaks$orange<-peaks$orange-correctOrange*peaks$green
      #   	}
      #   }
      
      
      ## MOVED TO qualityAssessment.R
      ### remove outliers using the control channel (e.g. blue)
      # if (!is.na(Cchannel)){
      #   lowerq = quantile(get(Cchannel, peaks))[2]
      #   upperq = quantile(get(Cchannel, peaks))[4]
      #   iqr = upperq - lowerq #Or use IQR(data)
      #   #mild.threshold.upper = (iqr * 1.5) + upperq
      #   #mild.threshold.lower = lowerq - (iqr * 1.5)
      #   extreme.threshold.upper = (iqr * 3) + upperq
      #   extreme.threshold.lower = lowerq - (iqr * 3)
      #   
      #   ixOut<-which((get(Cchannel, peaks))>extreme.threshold.upper | (get(Cchannel, peaks))<extreme.threshold.lower)
      #   if (length(ixOut)>0){peaks <- peaks[-ixOut,]}
      #   
      # }
    }
  }
  

  
  rownames(peaks)<-NULL
  return(peaks)
  
}
