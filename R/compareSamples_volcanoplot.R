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
#' Compare the different conditions across samples
#' 
#' \code{compareSamples_volcanoplot} 
#' 
#' This function compare the distribution of the values measures for each condition across 2 or more samples computing p-values
#' and effect size using parametric or non parametric tests. The comparison is visualised with a volcano plot showing the effect size vs the -log p-val
#' 
#' @param allData list with one element for each patient sample/cell line, there must be at least 2 runs (i.e. full cycle)
#' for each sample
#' @param controlName name of the control sample, default is "FS+FS"
#' @param thPval threshold on the p-value to be considered significant
#' @param thEffectSize threshold on the effect size
#' @importFrom lsr cohensD
#' @examples 
#' data("allData", package="BraDiPluS")
#' allData<-allData[c(1,2)]
#' res<-compareSamples(allData, controlName="FS + FS")
#' @export
#' 
compareSamples_volcanoplot <- function(allData, controlName="FS + FS", thPval=0.1, thEffectSize=0.3){
  
  library(ggplot2)
  library(ggrepel)
  
  # for each cell-line/patient
  res<-lapply(allData, function(myData){
    # compute median value for each run
    runs_medians<-do.call(cbind, lapply(myData, function(myRun){
      sapply(sapply(myRun, get, x="green"), median)
    }))
    # compute corresponding z-score
    runs_zscore<-apply(runs_medians,2,scale)
    rownames(runs_zscore)<-rownames(runs_medians)

    ## NOTE: for now I do it also on controls but I should consider removing them    
    # # find the control samples
    samplesNames<-rownames(runs_zscore)
    ix_control<- which(samplesNames==controlName)
    # controlData<-c(runs_zscore[ix_control,])
    # controlData<-controlData[!is.na(controlData)]
    # 
    # # remove the control samples
    runs_zscore<-runs_zscore[-ix_control,]
    # 

    ######
    is.control<- as.factor(rownames(runs_zscore)==controlName)
    nRuns<-ncol(runs_zscore)
    
    allData.df<-data.frame(ID=as.factor(rep(seq(1:nrow(runs_zscore)), nRuns)),
                           names=rep(rownames(runs_zscore),nRuns),
                           run=rep(colnames(runs_zscore), each=nrow(runs_zscore)),
                           values= do.call(rbind, lapply(split(runs_zscore, col(runs_zscore)), as.matrix)),
                           is.control=rep(is.control, nRuns))
    
    
    return(allData.df)
  })

  allData.df<-do.call(rbind, lapply(names(res), function(x){
    tmp<-res[[x]]
    tmp$patientCellLine=x
    return(tmp)
  }))

  res<-do.call(rbind, lapply(levels(allData.df$ID), function(x){
    tmp<-subset(allData.df, ID==x)
    
    ## tests for 2 groups (cell line/patient)
    parametric=T # for now only parametric test implemented: ToDo add non parametric and set this as argument with:
    #' TRUE to use parametric tests (t-test if 2 samples, ANOVA if more samples), FALSE to use non parametric
    #' #' tests (Wilcoxon rank sum test if 2 samples, Kruskal-Wallis rank sum test if more than 2 samples).
    if (length(allData)==2){
      if (parametric==T){  # t-test
        p.val = t.test(values ~ patientCellLine, data=tmp)$p.value
        effSize = cohensD(values ~ patientCellLine, data=tmp)
        # effSize = max(medianVals$values[2], medianVals$values[1])
        medianVals<-aggregate(values~patientCellLine, data=tmp, FUN=median)
        signR<-sign(medianVals$values[2] - medianVals$values[1])
      }else{  # Wilcoxon rank sum test
        # TODO
      }
      
    }else{ ## tests for more than 2 groups (cell line/patient)
      if (parametric==T){  # ANOVA
        # TODO
        
        signR=1 # here I don't care about sign
        
      }else{  # Kruskal-Wallis rank sum test
        # TODO
        
        signR=1 # here I don't care about sign
      }
    }
    
    
    return(data.frame(ID=x,
                      names=as.character(allData.df$names[which(allData.df$ID==x)[1]]),
                      p.val=p.val,
                      effSize=effSize,
                      signR=signR))
    
  }))
  
  # adjust p-values for multiple comparison testing
  res$p.val<-p.adjust(res$p.val, method = "BH")
  
  # define significance based on user defined thresholda
  res$threshold = as.numeric(as.factor(res$p.val <= thPval & abs(res$effSize) > thEffectSize))
  res$threshold<-factor(res$threshold, levels = c(1,2), labels = c("notSign", "Sign"))
  
  
  library(ggrepel)
  g = ggplot(data=res, aes(x=effSize*signR, y=-log10(p.val), colour=threshold, size=20), environment = environment()) +
    geom_point(alpha=0.5) +   theme(legend.position="none") +
    theme(legend.position = "none")+
    ylab("-log10 p-value") + ggtitle("") + xlab("z-score") + #xlab("cohen's D") +
    scale_colour_manual(values = c("notSign" = "#ededed", "Sign" = "#490A3D")) +
    theme_bw() + geom_hline(yintercept = -log10(thPval), linetype = "longdash", colour="#9e9e9e")+
    geom_vline(xintercept = thEffectSize, linetype = "longdash", colour="#9e9e9e")+
    geom_vline(xintercept = -thEffectSize, linetype = "longdash", colour="#9e9e9e")+
    geom_vline(xintercept = 0, colour="#9e9e9e") #+
  # xlim(-1, 1) #+ ylim(0, 95)  #xlim(-2, 2) + ylim(0, 85)
  g=g+geom_text_repel(data=subset(res, threshold!="notSign"), aes(x=effSize*signR, y=-log10(p.val), label=names), show.legend = NA, inherit.aes = F) +
    theme(legend.position = "none")
  plot(g)
  

  return(res)
}
  
  
  