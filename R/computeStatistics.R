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
#' Compute and plot z-score and p-value.
#' 
#' \code{computeStatistics} compute and plot z-score and p-value for one or more patient sample/cell lines.
#' 
#' This function computes the z-score and the p-values for all patient samples/cell lines provided as argument. For each 
#' patient sample/cell lines, the z-scores are computed separately for each run and then averaged. P-values are also computed
#' separately for each run, using Wilcoxon rank sum test to verify if the response to the tested drugs (alone or in combination)
#' is significantly higher than the measurements for control condition (no perturbation). P-values were FDR-corrected for multiple
#' hypotheses testing and combined across different runs using Fisher’s method.
#' 
#' @param allData list with one element for each patient sample/cell line, there must be at least 2 runs (i.e. full cycle)
#' for each sample
#' @param controlName name of the control sample, default is "FS+FS"
#' @param thPval threshold on the p-value to be considered significant
#' @param thSD threshold on the maximum allows standard deviation of z-scores for each sample across runs. Default is NA, so not sample
#' is discarted based on this; can be set, e.g. to 1, to discart samples with high variability across runs.
#' @param subsample select TRUE to consider always (for all patient sample/cell line) the same number of runs when computing the combined p-value,
#' this makes points in the volcano plot more comparable. FALSE all runs available are considered.
#' @param saveFiles TRUE to save the plot as pdf. FALSE to plot it.
#' @param showLabels TRUE to visualiza samples names.
#' @importFrom survcomp combine.test
#' @importFrom gridExtra grid.arrange
#' @examples 
#' data("allData", package="BraDiPluS")
#' allData<-list(patient3=allData[[5]])
#' res<-computeStatistics(allData, controlName="FS + FS", thPval=0.1, thSD=1, subsample=F, saveFiles=T)
#' @export
#' 
computeStatistics <- function(allData, controlName="FS + FS", thPval=0.1, thSD=NA, subsample=F, saveFiles=F, showLabels=F){
  
  # library(survcomp)
  # library(ggrepel)
  library(corrplot)
  library(ggplot2)
  library(gridExtra)
  library(corrplot) # note: this uses a slightly modified version of corrplot package (to handle NAs)
  if (showLabels==T){
    library(ggrepel)
  }
  
  # minimum number of runs (for subsampling when computing combined p-value)
  minRuns<-min(sapply(allData, length))
  
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
        if (nrow(x)>=3){ # needs at least 3 replicates
          wilcox.test(x$green, allControlRep, alternative="greater")$p.value
        }else{
          return(NA)
        }
      })
      thePval<-p.adjust(thePval, "BH")
      return(thePval)
    }))
    # combine p-values
    runs_pvalue_combined<-apply(runs_pvalue, 1, function(x){
      # NOTE: I could add as weight the sample size of each study
      # subsample runs when computing combined p-val (avoid diffeneces in volcano just due to number of runs)
      if (length(x)>minRuns & subsample==T){x<-sample(x,2, replace = F)}
      survcomp::combine.test(x, method = "fisher", na.rm = T)
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
    
    ######
    # prepare for boxplot
    is.significant<-as.factor(runs_zscore_median>0 & runs_pvalue_combined <= thPval)
    is.significant[ix_highSD]<-FALSE
    is.control<- as.factor(names(is.significant)==controlName)
    nRuns<-ncol(runs_zscore)
    dataBoxplot<-data.frame(ID=as.factor(rep(seq(1:nrow(runs_zscore)), nRuns)),
                            names=rep(names(is.significant),nRuns),
                            run=rep(colnames(runs_zscore), each=nrow(runs_zscore)),
                            values= do.call(rbind, lapply(split(runs_zscore, col(runs_zscore)), as.matrix)),
                            is.significant=rep(is.significant, nRuns),
                            is.control=rep(is.control, nRuns))
    
    # g_box <- ggplot(dataBoxplot, aes(x=ID, y=values, colour=is.significant, fill=is.control)) + geom_hline(yintercept=0, colour="grey70") +
    #   geom_boxplot(outlier.colour=NA) + 
    #   geom_point(position = position_jitter(width = 0.2)) + theme_bw() +
    #   scale_y_continuous(sec.axis = dup_axis(name = waiver())) +
    #   ylab("z-score") + xlab("") +
    #   scale_x_discrete(breaks = factor(seq(1:nrow(runs_zscore))), labels = names(is.significant)) +
    #   theme(axis.text.x=element_text(angle = -45, hjust = 0), legend.position = "top") +
    #   scale_colour_manual(values = c("grey80", "grey40"), name="", breaks=levels(is.significant), labels=c("not significant", "significant")) +
    #   scale_fill_brewer(palette="Accent", name="        ", breaks=levels(is.control), labels=c("drug response sample\n(single drug or pairwise combination)", "control sample\n(only FreeStyle medium, no drugs)"))

    ######
    # prepare for heatmap
    dataHeatmap<-data.frame(sample=samplesNames, zscore=runs_zscore_median, pvalue=runs_pvalue_combined)
    perturbations<-do.call(rbind, lapply(as.character(dataHeatmap$sample), function(theName){
      strsplit(theName, split = " \\+ ")[[1]]
    }))
    dataHeatmap<-cbind(perturbations, dataHeatmap, stringsAsFactors=FALSE)
    
    allCompounds<-unique(do.call(c, lapply(as.character(dataHeatmap$sample), function(theName){
      strsplit(theName, split = " \\+ ")[[1]]
    })))
    allCompounds[2:length(allCompounds)]<-sort(allCompounds[2:length(allCompounds)])
    
    dataHeatmap.m.zscore<-matrix(NA, nrow=length(allCompounds), ncol=length(allCompounds))
    colnames(dataHeatmap.m.zscore)<-rownames(dataHeatmap.m.zscore)<-allCompounds
    dataHeatmap.m.pval<-dataHeatmap.m.zscore
    
    for (i in 1:nrow(dataHeatmap)){
      # cat(dataHeatmap[i,1], dataHeatmap[i,2], "\n")
      dataHeatmap.m.zscore[dataHeatmap[i,1], dataHeatmap[i,2]]<-dataHeatmap.m.zscore[dataHeatmap[i,2], dataHeatmap[i,1]]<-dataHeatmap[i,"zscore"]
      dataHeatmap.m.pval[dataHeatmap[i,1], dataHeatmap[i,2]]<-dataHeatmap.m.pval[dataHeatmap[i,2], dataHeatmap[i,1]]<-dataHeatmap[i,"pvalue"]
    }
    
    is.control<- as.factor(dataHeatmap$sample==controlName)
    diag(dataHeatmap.m.zscore)<- -0.01
    dataHeatmap.m.zscore["FS", "FS"]<-median(dataHeatmap[is.control==T,"zscore"])
    # dataHeatmap.m.pval["FS", "FS"]<-0
    
    # to consider non significant also those with z-score<0
    diag(dataHeatmap.m.pval)<-0
    dataHeatmap.m.pval[dataHeatmap.m.zscore<=0]<-1
    
    ######
    # prepare for volcano plot
    dataVolcano<-data.frame(zscore=runs_zscore_median_rm, pvalue=runs_pvalue_combined_rm, sample=names(runs_pvalue_combined_rm))
    # dataVolcano$pvalue<-p.adjust(dataVolcano$pvalue, method = "BH")
    
    dataVolcano$threshold = as.numeric(as.factor(dataVolcano$pvalue <= thPval & dataVolcano$zscore >0))
    dataVolcano$threshold<-factor(dataVolcano$threshold, levels = c(1,2), labels = c("notSign", "Sign"))
    
    
    # g_vol = ggplot(data=dataVolcano, aes(x=zscore, y=-log10(pvalue), colour=threshold, size=10), environment = environment()) +
    #   geom_point(alpha=1) +   theme(legend.position="none") +
    #   theme(legend.position = "none") +
    #   # xlim(c(-2.8, 2.8)) + ylim(c(0, 2.6)) +
    #   xlab("effect size") + ylab("-log10 p-value") + ggtitle("") +
    #   scale_colour_manual(values = c("notSign" = "#ededed", "Sign" = "#74d600")) +
    #   theme_bw() + geom_hline(yintercept = -log10(0.1), linetype = "longdash", colour="#9e9e9e") +
    #   geom_vline(xintercept = 0, linetype = "longdash", colour="#9e9e9e")
    # 
    # if (showLabels==T){
    #   g_vol<-g_vol+geom_text_repel(data=subset(dataVolcano, threshold=="Sign"), aes(x=zscore, y=-log10(pvalue), label=sample), show.legend = NA, inherit.aes = F)
    # }
    # 
    # return(list(g_vol=g_vol, g_box=g_box, dataVolcano=dataVolcano, dataBoxplot=dataBoxplot))
    return(list(dataVolcano=dataVolcano, dataBoxplot=dataBoxplot, dataHeatmap.m.pval=dataHeatmap.m.pval, dataHeatmap.m.zscore=dataHeatmap.m.zscore))
  })
  

  ###
  # overall heatmap
  if (saveFiles==T){
    pdf("heatmap.pdf",width=4*length(allData), height=4,paper='special') 
    par(mfrow=c(1, length(res)))
  }
  
  cl.lim<-c(-2,2.5)
  cols<-colorRampPalette(c("grey60", "grey70", "grey80", '#FFFFFF', '#DB949D', "#B61126", '#b20319'))(20)
  for (x in names(res)){
    tmp.zscore<-res[[x]]$dataHeatmap.m.zscore
    tmp.zscore[tmp.zscore<cl.lim[1]]<-cl.lim[1]
    tmp.zscore[tmp.zscore>cl.lim[2]]<-cl.lim[2]
    corrplot(tmp.zscore, insig="pch", pch="x", pch.cex=1, is.corr = F,  method="color", na.label = "o", cl.length=10, title = x,
             p.mat = res[[x]]$dataHeatmap.m.pval, sig.level=thPval, col=cols, cl.lim=c(-2,2.5), tl.col="black", addgrid.col="grey80", cl.align.text="l")
  }
  
  if (saveFiles==T){
    dev.off()
  }else{
    readline("press Enter to continue...")
  }
  
  
  
  
  
  
  
  ###
  # overall boxplot
  
  allBoxplot.df<-do.call(rbind, lapply(names(res), function(x){
    tmp<-res[[x]]$dataBoxplot
    tmp$patientCellLine=x
    return(tmp)
  }))

  gBoxAll <- ggplot(allBoxplot.df, aes(x=ID, y=values, colour=is.significant, fill=is.control)) + geom_hline(yintercept=0, colour="grey70") +
    geom_boxplot(outlier.colour=NA) + 
    geom_point(position = position_jitter(width = 0.2)) + theme_bw() +
    scale_y_continuous(sec.axis = dup_axis(name = waiver())) +
    ylab("z-score") + xlab("") +
    scale_x_discrete(breaks = as.numeric(levels(allBoxplot.df$ID)), labels = allBoxplot.df$names[as.numeric(levels(allBoxplot.df$ID))]) +
    theme(axis.text.x=element_text(angle = -45, hjust = 0), legend.position = "top") +
    scale_colour_manual(values = c("grey80", "grey40"), name="", breaks=levels(allBoxplot.df$is.significant), labels=c("not significant", "significant")) +
    scale_fill_brewer(palette="Accent", name="        ", breaks=levels(allBoxplot.df$is.control), labels=c("drug response sample\n(single drug or pairwise combination)", "control sample\n(only FreeStyle medium, no drugs)")) +
    facet_grid(patientCellLine ~ .)
  if (saveFiles==T){
    ggsave(gBoxAll, file="boxplot.pdf", width = 12, height = 4*length(res))
  }else{
    print(gBoxAll)
    readline("press Enter to continue...")
  }
  
  
  
  ###
  # overall volcano plot
  
  allVolplot.df<-do.call(rbind, lapply(names(res), function(x){
    tmp<-res[[x]]$dataVolcano
    tmp$patientCellLine=x
    return(tmp)
  }))
  
  gVolAll = ggplot(data=allVolplot.df, aes(x=zscore, y=-log10(pvalue), colour=patientCellLine), size=5, environment = environment()) +
    geom_point(alpha=1) +
    xlab("z-score") + ylab("-log10 p-value") + ggtitle("") +
    theme_bw() + geom_hline(yintercept = -log10(0.1), linetype = "longdash", colour="#9e9e9e") +
    geom_vline(xintercept = 0, linetype = "longdash", colour="#9e9e9e") +
    scale_colour_brewer(palette="Set2", name="")

  if (showLabels==T){
    gVolAll<-gVolAll+geom_text_repel(data=subset(allVolplot.df, threshold=="Sign"), aes(x=zscore, y=-log10(pvalue), label=sample), size = 3, col="grey60", show.legend = NA, inherit.aes = F) 
  }
  
  if (saveFiles==T){
    ggsave(gVolAll, file="volcanoplot.pdf", width = 6, height = 5)
  }else{
    gVolAll
  }
  

  
  return(allVolplot.df)
  # return(list(myVolcanoPlot=gVolAll, myBoxPlot=gBoxAll, dataVolcano=allVolplot.df, dataBoxplot=allBoxplot.df))
}
  
  
  