library(BraDiPluS)
# data("allData", package="BraDiPluS")
# load("/Users/eduati/BraDiPluS/data/allData_v2.RData")

# load("/Users/eduati/BraDiPluS/data/allData.RData")
# discarded<-allData[c(6)]
# rm(allData)
load("/Users/eduati/BraDiPluS/data/allData_v2.RData")

myData<-allData[[2]]

res_all<-lapply(1:length(myData), function(myRun_ix){
  myRun<-myData[[myRun_ix]]
  res<-do.call(rbind, lapply(1:length(myRun), function(myCondition_ix){
    # cat(myCondition_ix, "\n")
    myCondition<-myRun[[myCondition_ix]]
    if (nrow(myCondition)>0){
      res_tmp<-data.frame(blue_value=myCondition$blue, sample=myCondition_ix, count=seq(1:nrow(myCondition)))
    }else{
      res_tmp<-data.frame(blue_value=NULL, sample=NULL, count=NULL)
    }
    return(res_tmp)
  }))
  res<-res[-1,]
  res$run<-myRun_ix
  return(res)
})

# x<-do.call(rbind, res_all)

my_summary<-do.call(cbind, lapply(res_all, function(x){
  
  # lower_outlier<-quantile(x$blue_value)[2]-1.5*IQR(x$blue_value)
  # upper_outlier<-quantile(x$blue_value)[4]+1.5*IQR(x$blue_value)
  upper_outlier<-quantile(x$blue_value, probs=c(0.95))
  
  tmp<-sapply(sort(unique(x$count)), function(y){
    tmp<-subset(x, count==y)
    sum(tmp$blue_value>upper_outlier)*100/nrow(tmp)
  })
  
  tmp<-tmp[1:8]
  
  return(tmp)
  # do.call(c, lapply(sort(unique(x$count)), function(y){
  #   wilcox.test(subset(x, count==y)$blue_value, res_all$blue_value, alternative = "greater")$p.value
  # }))
}))

apply(my_summary,1,median)

# res_all$label<-(res_all$count==1)
# # res_all<-lapply(res_all, function(x){x$label <- x$count==1; return(x)})
# 
# library(ggplot2)
# ggplot(res_all[[1]], aes(x=blue_value, fill=label)) + geom_density(alpha=.3)
# ggplot(res_all[[1]], aes(x=blue_value, fill=as.factor(count))) + geom_density(alpha=.3)
# 
# wilcox.test(subset(res_all[[1]], count==9)$blue_value, res_all$blue_value, alternative = "greater")
# wilcox.test(subset(res_all, count==5)$blue_value, res_all$blue_value, alternative = "greater")
# 
# 
# ggplot(res_all, aes(x=count, y=blue_value)) + 
#   geom_boxplot() + facet_grid(run ~ ., scales="free") + theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(x="", y = "")  + geom_jitter(width = 0.2)
# 
# 
# ggplot(res_all[[1]], aes(x=blue_value)) + geom_density(alpha=.3)

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
# patients<-c(discarded, patients)
# names(patients)[1]<-"biopsy: discarded sample"

# compute and plot z-score for all cell line and patient samples 
res_stats<-computeStatistics(allData = patients, controlName="FS + FS", thPval=0.05, thSD=NA, subsample=F, saveFiles=F)

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
write.csv(mySpiderData, file="spiderData.csv", quote = F)


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


allConditions<-sort(levels(res_stats$sample))
allSamples<-levels(res_stats$patientCellLine)

matrixSign<-matrix(NA, ncol=length(allConditions), nrow=length(allSamples))
colnames(matrixSign)<-allConditions
rownames(matrixSign)<-allSamples
matrixValue<-matrix(NA, ncol=length(allConditions), nrow=length(allSamples))
colnames(matrixValue)<-allConditions
rownames(matrixValue)<-allSamples

for (i in 1:nrow(res_stats)){
  if (is.na(res_stats$threshold[i])){
    
  }
  else if (as.character(res_stats$threshold[i])=="Sign"){
    matrixSign[res_stats$patientCellLine[i], res_stats$sample[i]]<-1
  }else if (as.character(res_stats$threshold[i])=="notSign"){
    matrixSign[res_stats$patientCellLine[i], res_stats$sample[i]]<-0
  }
  
  matrixValue[res_stats$patientCellLine[i], res_stats$sample[i]]<-res_stats$zscore[i]
}



all_drugs<-unique(do.call(c, lapply(allConditions, function(x){
  tmp<-strsplit(x, " \\+ ")[[1]]
  tmp<-sub(" ", "", tmp)
  return(tmp)
})))
all_drugs<-sort(setdiff(all_drugs, "FS"))

all_targets<-list()
all_targets[[1]]<-"IKK"
all_targets[[2]]<-"MEK"
all_targets[[3]]<-"JAK"
all_targets[[4]]<-"PI3K"
all_targets[[5]]<-"EGFR"
all_targets[[6]]<-"DNA replication"
all_targets[[7]]<-"AKT"
all_targets[[8]]<-"DNA replication"
all_targets[[9]]<-c("AKT", "PDPK1")
all_targets[[10]]<-c("TNF")

names(all_targets)<-all_drugs


allTargets<-sort(unique(unlist(all_targets)))
allTargets<-allTargets[c(2, 1, 3:9)]
matrixTargets<-matrix(0, ncol=length(allConditions), nrow=length(allTargets))
colnames(matrixTargets)<-allConditions
rownames(matrixTargets)<-allTargets


for (x in allConditions){
  tmp<-strsplit(x, " \\+ ")[[1]]
  tmp<-sub(" ", "", tmp)
  zz<-unname(unlist(all_targets[tmp]))
  matrixTargets[zz, x]<-1
}



ix_order<-order(apply(matrixSign,2,function(x){sum(x, na.rm = T)}), decreasing = T)
corrplot(matrixSign[,ix_order], is.corr = F, na.label = "o", tl.col="black", addgrid.col="white", cl.align.text="l", pch.cex=1) 

corrplot(matrixTargets[,ix_order], is.corr = F, na.label = "o", tl.col="black", addgrid.col="grey80", cl.align.text="l", pch.cex=1, method="square") 


matrixValue[,ix_order]

## discarded patient
patients<-allData[c(6)]

# compute and plot z-score for all cell line and patient samples 
res_stats<-computeStatistics(allData = patients, controlName="FS + FS", thPval=0.05, thSD=NA, subsample=F, saveFiles=T)
