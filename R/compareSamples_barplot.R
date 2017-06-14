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
#' \code{compareSamples_barplot} 
#' 
#' This function compare the distribution of the values measures for each condition across 2 or more samples. 
#' For each cell line/patient it tests if each experimental condition is more promising than in the compared cell line/patient
#' and visualises with a barplot the conditions which are significantly different ordered from most to least promising (based on z-score)
#' 
#' @param allData list with one element for each patient sample/cell line, there must be at least 2 runs (i.e. full cycle)
#' for each sample
#' @param controlName name of the control sample, default is "FS+FS"
#' @param thPval threshold on the p-value to be considered significant
#' @param res_stats output of function computeStatistics, to consider only significant conditions. Default is NA.
#' @importFrom lsr cohensD
#' @examples 
#' data("allData", package="BraDiPluS")
#' allData<-allData[c(1,2)]
#' res<-compareSamples_barplot(allData, controlName="FS + FS")
#' @export
#' 
compareSamples_barplot <- function(allData, controlName="FS + FS", thPval=0.05, res_stats=NULL){
  
  library(ggplot2)
  

  # for each cell-line/patient
  res<-lapply(allData, function(myData){
    # compute median value for each run
    runs_medians<-do.call(cbind, lapply(myData, function(myRun){
      sapply(sapply(myRun, get, x="green"), median)
    }))
    # compute corresponding z-score
    runs_zscore<-apply(runs_medians,2,scale)
    rownames(runs_zscore)<-rownames(runs_medians)

    
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

  # remove the control samples
  allData.df<-subset(allData.df, is.control==F)
  
  res<-do.call(rbind, lapply(levels(droplevels(allData.df$ID)), function(x){
    
    tmp<-subset(allData.df, ID==x)
    medianVals<-aggregate(values~patientCellLine, data=tmp, FUN=median, na.action = na.omit)

    res2<-do.call(rbind, lapply(medianVals$patientCellLine, function(y){
      refGroup<-subset(tmp, patientCellLine==y)
      medianValue=median(refGroup$values, na.rm = T)
      # hist(refGroup$values, main=y)
      res3<-do.call(rbind, lapply(setdiff(medianVals$patientCellLine, y), function(yy){
        compGroup<-subset(tmp, patientCellLine==yy)
        # hist(compGroup$values, main=yy)
        
        p.val=wilcox.test(refGroup$values, compGroup$values, alternative="greater")$p.value
        # p.val=t.test(refGroup$values, compGroup$values, alternative="greater")$p.value
        
        effSize=10   ## fon now I am not considering effectSize: put here just arbitrarily high number
        
        # ## alternative with effectSize
        # library(coin)
        # my.df<-rbind(refGroup, compGroup)
        # my.df$patientCellLine<-factor(my.df$patientCellLine)
        # wt_res<-wilcox_test(values ~ patientCellLine, data = my.df, alternative = "greater", distribution="exact")
        # effSize<-wt_res@statistic@teststatistic/sqrt(nrow(my.df))
        # wt_res@statistic
        # 
        ## otherwise with parametric test
        # p.val=t.test(refGroup$values, compGroup$values, alternative="greater")$p.value
        # effSize=cohensD(refGroup$values, compGroup$values)
        return(data.frame(effSize=effSize, p.val=p.val))
      }))
      thEffectSize<-0 ## fon now I am not considering effectSize: put here 0
      is.signif.comparison=any(res3$effSize>thEffectSize & res3$p.val<thPval)
      if (!is.null(res_stats)){  #if results from computeStatistic are availe I mark significant conditions 
        myPatinetCellLine<-which(res_stats$patientCellLine==y)
        myCondition<-which(res_stats$sample==as.character(unique(tmp$names)))
        ix<-intersect(myCondition, myPatinetCellLine)
        is.signif.condition=res_stats[ix,"threshold"]
      }else{
        is.signif.condition=NA
      }
      
      return(data.frame(medianValue=medianValue, is.signif.comparison=is.signif.comparison, is.signif.condition=is.signif.condition, patientCellLine=y))
    }))
    
    res2$ID=x
    res2$names=unique(tmp$names)
    
    return(res2)
  }))
  
  res_all<-res
  
  # if I want to plot only the significant ones
  res<-subset(res, is.signif.comparison==T) # comment this to show also non significant comparisons
  if (!is.null(res_stats)){ # if info available consider only significant conditions for each cell line/patient
    res<-subset(res, is.signif.condition=="Sign")
  }
  
  g<-lapply(levels(res$patientCellLine), function(x){
    tmp<-subset(res, patientCellLine==x)
    tmp<-tmp[order(tmp$medianValue),]
   
    ## this is to be used if I want to plot also non significant ones
    # gg <- ggplot(tmp, aes(x=reorder(names, medianValue), y=medianValue)) + 
    #   geom_point(size=3, aes(colour=as.factor(is.significant))) + 
    #   geom_segment(aes(xend=reorder(names, medianValue)), yend=0) +
    #   ylab("z-score") + xlab("conditions") + ggtitle("") +
    #   scale_colour_manual(values = c("grey60", '#b20319')) +
    #   expand_limits(y=0) +
    #   coord_flip() + theme_bw() + geom_hline(yintercept = 0, linetype = "longdash", colour="#9e9e9e")
    # 
    gg <- ggplot(tmp, aes(x=reorder(names, medianValue), y=medianValue)) + 
      geom_point(size=3, col='#b20319') + 
      geom_segment(aes(xend=reorder(names, medianValue)), yend=0) +
      ylab("z-score") + xlab("promising\nconditions") + ggtitle(x) +
      expand_limits(y=0) +
      coord_flip() + theme_bw() + geom_hline(yintercept = 0, linetype = "longdash", colour="#9e9e9e")

    ## this to show different colour for significant and not (only if I am plotting also non significant comparisons)   
    # gg <- ggplot(tmp, aes(x=reorder(names, medianValue), y=medianValue, col=is.signif.comparison)) + 
    #   geom_point(size=3) + 
    #   geom_segment(aes(xend=reorder(names, medianValue)), yend=0) +
    #   ylab("z-score") + xlab("promising\nconditions") + ggtitle(x) +
    #   expand_limits(y=0) +
    #   coord_flip() + theme_bw() + geom_hline(yintercept = 0, linetype = "longdash", colour="#9e9e9e")
    # 
    
    return(gg)  
  })
  
  library(gridExtra)
  do.call("grid.arrange", c(g, ncol=1))
  
  return(res_all)
}
  
  
  