# INPUT:
# runs: list with each run as one element
# OUTPUT:
# runs: same list but replicates that are outliers based on orange channel are removed, and samples with less then 2 replicates are also removed

# run for each run the removeOutliers
qualityAssessment <- function(runs){
  par(mfrow=c(1,length(runs)))
  
  res<-list()
  for (i in 1:length(runs)){
    res[[i]]<-removeOutliers(runs[[i]], i)
  }
  
  return(res)

}


# remove ouliers for each run
removeOutliers <- function(x, index){
  samples.names<-names(x)
  
  allData<-unlist(sapply(x, get, x="orange"))

  bxp<-boxplot(allData, ylab="orange (all replicates, all samples)")
  main<-paste("run", index, ": median=", round(median(allData), 2), ", IQR=", round(IQR(allData), 3))
#   title(main=main, sub=paste("total =", length(allData), "   outliers =", round(length(bxp$out),3)))
  
  lower_outlier<-quantile(allData)[2]-1.5*IQR(allData)
  upper_outlier<-quantile(allData)[4]+1.5*IQR(allData)
  
  # remove outliers from each sample
  x_new<-lapply(x, function(y){
    subset(y, (orange > lower_outlier & orange < upper_outlier))
  })

  # find samples with less then 2 replicates after removing outliers
  nRepl<-sapply(x_new, nrow)
  ixRemove<-which(nRepl<2)
  
  cat("\nrun", index, ":\n")
  cat("Sample(s): ", paste(names(x_new)[ixRemove], sep=","), "have less then 2 replicates and will therefore be removed")

  if (length(ixRemove)>0){
    sub1<-paste("samples removed:", paste(names(x_new)[ixRemove], collapse=","))
  }else{
    sub1<-"samples removed: none"
  }
  
  sub2<-paste("total =", length(allData), "   outliers =", round(length(bxp$out),3))
  title(main=main, sub=paste(sub1, sub2, sep="\n"))
  
  
#   if (length(ixRemove)>0){
#     x_new<-x_new[-ixRemove]
#   }
  
  # instead of removing it, I substitute it with an empty data frame
  for (i in ixRemove){
    x_new[[i]]<-setNames(data.frame(matrix(ncol = ncol(x_new[[i]]), nrow = 0)), colnames(x_new[[i]]))
  }

  return(x_new)
}


