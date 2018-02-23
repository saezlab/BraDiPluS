# #########################################################################################################
# Script to analyis data from mouse model experiments
# #########################################################################################################

# ****************
# packages
library(ggplot2)
library(gdata)
library(reshape2)

# ****************
# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ****************
# working directory
AsPC1 <- read.xls('../data/mouse_model_22feb18.xls', sheet = 1, stringsAsFactors=FALSE)
AsPC1.df <- AsPC1[4:nrow(AsPC1), 3:ncol(AsPC1)]
for (i in 1:nrow(AsPC1.df)){
  AsPC1.df[i,] <- (as.numeric(AsPC1.df[i,])-as.numeric(AsPC1.df[i,1]))/as.numeric(AsPC1.df[i,1])*100
}
rownames(AsPC1.df) <- AsPC1[4:nrow(AsPC1), 2]
colnames(AsPC1.df) <- AsPC1[3, 3:ncol(AsPC1)]
AsPC1.df$id <- rownames(AsPC1.df)
AsPC1.df <- melt(AsPC1.df, id.vars = "id")

BxPC3 <- read.xls('../data/mouse_model_22feb18.xls', sheet = 2, stringsAsFactors=FALSE)
BxPC3.df <- BxPC3[4:nrow(BxPC3), 3:ncol(BxPC3)]
rownames(BxPC3.df) <- BxPC3[4:nrow(BxPC3), 2]
colnames(BxPC3.df) <- BxPC3[3, 3:ncol(BxPC3)]
BxPC3.df$id <- rownames(BxPC3.df)
BxPC3.df <- melt(BxPC3.df, id.vars = "id")



