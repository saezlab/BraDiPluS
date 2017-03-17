data("allData", package="BraDiPluS")
load("/Users/eduati/BraDiPluS/data/allData.RData")
setwd("/Users/eduati/BraDiPluS/fig/")

# compute median value across replicate plugs for each run
allData_median<-lapply(allData, function(x){
  runs_medians<-do.call(cbind, lapply(x, function(myRun){
    sapply(sapply(myRun, get, x="green"), median)
  }))
})

# compute z-score
allData_scale<-lapply(allData_median, function(x){
  tmp<-apply(x,2,scale)
  rownames(tmp)<-rownames(x)
  return(tmp)
})
  
# boxplot with dots
# pdf("boxplots_runs.pdf",width=6,height=20,paper='special') 
par(mfrow=c(length(allData_scale),1),mar=c(15,5,3,3), cex.axis=0.7)
aa<-lapply(names(allData_scale), function(x){
  data<-allData_scale[[x]]
  data_scale.l <- setNames(split(data, seq(nrow(data))), rownames(data))
  col1<-rep("#b2d2e2", length(data_scale.l))
  col1[which(names(data_scale.l)=="FS + FS")]<-"#9998bb"
  col2<-rep("#4491b6", length(data_scale.l))
  col2[which(names(data_scale.l)=="FS + FS")]<-"#4c4a6f"
  boxplot(data_scale.l, col=col1, las=2, outline=F, main=x, xlab="", ylab="Casp3 activity")
  stripchart(data_scale.l, vertical=T, 
             method = "jitter", add = TRUE, pch = 20, col = col2)
  axis(side = 4)
})
# dev.off()


load("/Users/eduati/BraDiPluS/data/allData.RData")
source('~/BraDiPluS/R/volcanoPlotResponsiveSamples.R')
myData<-allData[[1]]
thSD<-1.75
library(survcomp)
volcanoPlotResponsiveSamples(allData, controlName="FS + FS", thSD=NA)
  



