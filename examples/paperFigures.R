library(BraDiPluS)
# data("allData", package="BraDiPluS")
load("/Users/eduati/BraDiPluS/data/allData_v2.RData")
setwd("/Users/eduati/BraDiPluS/fig/")

res_stats<-computeStatistics(allData = allData, controlName="FS + FS", thPval=0.05, thSD=NA, subsample=F, saveFiles=F)


###
## only for cell lines
cellLines<-allData[c(1,2)]

# compute and plot z-score for all cell line and patient samples 
res_stats<-computeStatistics(allData = cellLines, controlName="FS + FS", thPval=0.05, thSD=NA, subsample=F, saveFiles=F)


# compare conditions between samples 
# res<-compareSamples_volcanoplot(allData = cellLines, controlName="FS + FS", thPval=0.1, thEffectSize=0.3)
res<-compareSamples_barplot(allData = cellLines, controlName="FS + FS", thPval=0.05, res_stats=res_stats)


###
## only for patients
# patients<-allData[c(3,4,5,6)]
patients<-allData[c(4,6,5,3)]

# compute and plot z-score for all cell line and patient samples 
res_stats<-computeStatistics(allData = patients, controlName="FS + FS", thPval=0.05, thSD=NA, subsample=F, saveFiles=T)

# compare conditions between samples 
# res<-compareSamples_volcanoplot(allData = cellLines, controlName="FS + FS", thPval=0.1, thEffectSize=0.3)
res<-compareSamples_barplot(allData = patients, controlName="FS + FS", thPval=0.05, res_stats=res_stats)


res_stats$patientCellLine<-factor(res_stats$patientCellLine, levels = names(patients), labels = names(patients) , ordered=T)
g<-lapply(levels(res_stats$patientCellLine), function(x){
  tmp<-subset(res_stats, (patientCellLine==x & threshold=="Sign"))
  tmp<-tmp[order(tmp$zscore, decreasing = T),]
  
  tmp$otherNonSignSamples<-as.factor(sapply(tmp$sample, function(xx){
    tmp2<-subset(res_stats, sample==xx)
    return(sum(tmp2$threshold=="notSign", na.rm = T))
  }))
  levels(tmp$otherNonSignSamples)<-c(0,1,2,3)
  
  tmp<-tmp[1:10,]
  
  print(tmp)
  
  # this to show different colour for significant and not (only if I am plotting also non significant comparisons)
  gg <- ggplot(tmp, aes(x=reorder(sample, zscore), y=zscore, col=otherNonSignSamples)) +
    geom_point(size=3) +
    geom_segment(aes(xend=reorder(sample, zscore)), yend=0) +
    ylab("z-score") + xlab("promising\nconditions") + ggtitle(x) +
    expand_limits(y=0) +
    coord_flip() + theme_bw() + geom_hline(yintercept = 0, linetype = "longdash", colour="#9e9e9e") +
    # scale_color_brewer(palette="Blues") +
    scale_color_manual(limits = c(0,1,2,3),
                        values=c("#EEEEEE", "#BBBBBB", "#777777", "#222222")) + theme(legend.position="bottom") 
    
  
  
  return(gg)  
})


library(gridExtra)
do.call("grid.arrange", c(g, ncol=length(g)))


### try spider plot by drug
allCompounds<-sort(unique(do.call(c, lapply(levels(res_stats$sample), function(x){
  strsplit(x, split = " \\+ ")[[1]]
}))))
allCompounds<-setdiff(allCompounds, "FS")

res_stats$patientCellLine<-as.factor(res_stats$patientCellLine)

mySpiderData<-do.call(rbind, lapply(levels(res_stats$patientCellLine), function(x){
  tmp<-subset(res_stats, (patientCellLine==x))
  
  res<-sapply(allCompounds, function(y){
    median(tmp[grep(y, tmp$sample),"zscore"], na.rm = T)
  })
  
  return(res)
}))
rownames(mySpiderData)<-levels(res_stats$patientCellLine)
mySpiderData<-as.data.frame(mySpiderData)

library(fmsb)
generalMin<-floor((min(mySpiderData, na.rm = T))*10)/10
generalMax<-ceiling((max(mySpiderData, na.rm = T))*10)/10
mySpiderData=rbind(rep(generalMax,ncol(mySpiderData)) , rep(generalMin,ncol(mySpiderData)) , mySpiderData)

# radarchart(mySpiderData)

colors_border=c( rgb(108,167,93, 230, maxColorValue = 255),
                 rgb(148,115,198, 230, maxColorValue = 255),
                 rgb(191,130,59, 230, maxColorValue = 255),
                 rgb(204,84,109, 230, maxColorValue = 255) )
colors_in=c( rgb(108,167,93, 100, maxColorValue = 255),
             rgb(148,115,198, 100, maxColorValue = 255),
             rgb(191,130,59, 100, maxColorValue = 255),
             rgb(204,84,109, 100, maxColorValue = 255) )

radarchart( mySpiderData  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(mySpiderData[2,1],mySpiderData[1,1],length.out = 5), cglwd=0.8,
            #custom labels
            vlcex=0.8 
)
legend(x=0.9, y=1.35, legend = rownames(mySpiderData[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)

