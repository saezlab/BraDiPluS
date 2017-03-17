
volcanoPlotResponsiveSamples <- function(allData, controlName="FS + FS", thSD=NA){
  # for each cell-line/patient
  res<-lapply(allData, function(myData){
    # compute median value for each run
    runs_medians<-do.call(cbind, lapply(myData, function(myRun){
      sapply(sapply(myRun, get, x="green"), median)
    }))
    # compute p-value for each run (compared with control distribution)
    runs_pvalue<-do.call(cbind, lapply(myData, function(myRun){
      ix_control<- which(names(myRun)==controlName)
      allControlRep<-do.call(c, lapply(myRun[ix_control], function(x){
        x$green
      }))
      thePval<-sapply(myRun, function(x){
        if (nrow(x)>=3){
          wilcox.test(x$green, allControlRep, alternative="greater")$p.value
        }else{
          return(NA)
        }
      })
    }))
    # combine p-values
    runs_pvalue_combined<-apply(runs_pvalue, 1, function(x){
      # NOTE: I could add as weight the sample size of each study
      if (length(x)>2){x<-sample(x,2, replace = F)}
      combine.test(x, method = "fisher", na.rm = T)
      # mean(x, na.rm = T)
    })
    # compute corresponding z-score
    runs_zscore<-apply(runs_medians,2,scale)
    rownames(runs_zscore)<-rownames(runs_medians)
    # and the median z-score
    runs_zscore_median<-apply(runs_zscore, 1, function(x) median(x, na.rm = T))
    runs_zscore_sd<-apply(runs_zscore, 1, function(x) sd(x, na.rm = T))
    
    # find the control samples
    samplesNames<-rownames(runs_zscore)
    ix_control<- which(samplesNames==controlName)
    controlData<-c(runs_zscore[ix_control,])
    controlData<-controlData[!is.na(controlData)]
    
    # and remove those with with high variability among replicates I don't consider the data (put NA)
    if (!is.na(thSD)){
      ix_highSD<-which(runs_zscore_sd>thSD)
    }else{
      ix_highSD<-c()
    }
    
    # remove the control and the samples with high SD
    ix_remove<-c(ix_control, ix_highSD)
    runs_zscore_median_rm<-runs_zscore_median[-ix_remove]
    runs_pvalue_combined_rm<-runs_pvalue_combined[-ix_remove]
    dataVolcano<-data.frame(zscore=runs_zscore_median_rm, pvalue=runs_pvalue_combined_rm, sample=names(runs_pvalue_combined_rm))
    dataVolcano$pvalue<-p.adjust(dataVolcano$pvalue, method = "BH")
    
    dataVolcano$threshold = as.numeric(as.factor(dataVolcano$pvalue <= 0.1 & dataVolcano$zscore >0))
    dataVolcano$threshold<-factor(dataVolcano$threshold, levels = c(1,2), labels = c("notSign", "Sign"))
    
    library(ggrepel)
    g = ggplot(data=dataVolcano, aes(x=zscore, y=-log10(pvalue), colour=threshold, size=10), environment = environment()) +
      geom_point(alpha=1) +   theme(legend.position="none") +
      theme(legend.position = "none") +
      # xlim(c(-2.8, 2.8)) + ylim(c(0, 2.6)) +
      xlab("effect size") + ylab("-log10 p-value") + ggtitle("") +
      scale_colour_manual(values = c("notSign" = "#ededed", "Sign" = "#74d600")) +
      theme_bw() + geom_hline(yintercept = -log10(0.1), linetype = "longdash", colour="#9e9e9e") +
      geom_vline(xintercept = 0, linetype = "longdash", colour="#9e9e9e")
    
    return(list(g=g, dataVolcano=dataVolcano))
  })
  
  allRes.df<-do.call(rbind, lapply(names(res), function(x){
    tmp<-res[[x]]$dataVolcano
    tmp$patientCellLine=x
    return(tmp)
  }))
  
  library(ggrepel)
  # g = ggplot(data=allRes.df, aes(x=zscore, y=-log10(pvalue), colour=threshold, group=patientCellLine, size=10), environment = environment()) +
  #   geom_point(alpha=1) +   theme(legend.position="none") +
  #   theme(legend.position = "none") +
  #   # xlim(c(-2.8, 2.8)) + ylim(c(0, 2.6)) +
  #   xlab("effect size") + ylab("-log10 p-value") + ggtitle("") +
  #   scale_colour_manual(values = c("notSign" = "#ededed", "Sign" = "#74d600")) +
  #   theme_bw() + geom_hline(yintercept = -log10(0.1), linetype = "longdash", colour="#9e9e9e") +
  #   geom_vline(xintercept = 0, linetype = "longdash", colour="#9e9e9e")
  # 
  g = ggplot(data=allRes.df, aes(x=zscore, y=-log10(pvalue), colour=patientCellLine, size=10), environment = environment()) +
    geom_point(alpha=1) +   theme(legend.position="none") +
    theme(legend.position = "none") +
    # xlim(c(-2.8, 2.8)) + ylim(c(0, 2.6)) +
    xlab("z-score") + ylab("-log10 p-value") + ggtitle("") +
    theme_bw() + geom_hline(yintercept = -log10(0.1), linetype = "longdash", colour="#9e9e9e") +
    geom_vline(xintercept = 0, linetype = "longdash", colour="#9e9e9e")
  g
  return(g)
  
}
  
  
  