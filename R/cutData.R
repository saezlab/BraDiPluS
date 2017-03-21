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
#' Cut a portion of the data.
#' 
#' \code{cutData} returns a data structure limited to the selected time range.
#' 
#' This function allows to create a new data structure including only data in a user
#' defined time interval (e.g. to remove test measuremets or to select one
#' specific cycle).
#' 
#' @param data as imported from the .txt file saved by the in-house LabVIEW progam
#' @param startTime initial time point to be considered
#' @param endTime final time point to be considered
#' @return A data frame with the same structure of the input data but limited to
#' the selected time range.
#' @examples 
#' data(BxPC3_data, package="BraDiPluS")
#' MyDataNew <- cutData(data = MyData, startTime=6000, endTime=6500)
#' @export

cutData <- function(data, startTime=NA, endTime=NA){
  if (is.na(startTime)){startTime<-data$time[1]}
  if (is.na(endTime)){endTime<-data$time[nrow(data)]}
  
  ix1 <- which(data$time>=startTime)
  ix2 <- which(data$time>=endTime)
  
  data <- data[ix1[1]:ix2[1],]
  return(data)
}