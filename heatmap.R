###### filtering the list showing the heatmap (original code is from ‘NPY’ project)
rm(list=ls())
library(dplyr)
options(stringsAsFactors = FALSE)
### the function used to draw the heatmap
heatmapTCF <- function(inputA, inputB, figfile, mnames = c("DataA", "DataB"), minfold=0.8, maxfold=1.25, orderA = TRUE) {
  library(ggplot2)
  library(cowplot)
  library(reshape2)
  ## set column names for the inputA and inputB
  colnames(inputA) <- c("CycID", "Pva", "Pha", paste("sample", 1:(ncol(inputA) -3), sep=""))
  colnames(inputB) <- c("CycID", "Pva", "Pha", paste("sample", 1:(ncol(inputB) -3), sep=""))
  ## the sub function used to prepare the data frame for drawing heatmap
  figdataF <- function(inputD, minfold = minfold, maxfold = maxfold)  {
    ## arrange the input data by phase information
    inputD <- arrange(inputD, Pha)
    ## median normalization of the expression profiles
    figM <- as.matrix(inputD[, -(1:3)])
    fig_median <- apply(figM, 1, median)
    figM <- figM / fig_median
    figM[figM < minfold] <- minfold
    figM[figM > maxfold] <- maxfold
    figD <- as.data.frame(figM)
    ## prepare the data frame for heatmap
    id_order <- as.character(inputD$CycID)
    id_factor <- as.factor(length(id_order):1)
    figD <- mutate(figD, id_factor = id_factor)
    cnum <- ncol(figD)
    figD <- figD[,c(cnum, 1:(cnum - 1))]
    figD.m <- melt(figD)
    return(figD.m)
  }
  gid <- inputA$CycID
  rownames(inputA) <- inputA$CycID
  rownames(inputB) <- inputB$CycID
  ## put the same gene order between inputA and inputB
  inputB <- inputB[gid,]
  ## ignore the phase order of inputB, following the same phase order as inputA
  if (orderA) {
    inputB$Pha <- inputA$Pha
  }  else  {
    inputA$Pha <- inputB$Pha
  }
  figA.m <- figdataF(inputD = inputA, minfold = minfold, maxfold = maxfold)
  figB.m <- figdataF(inputD = inputB, minfold = minfold, maxfold = maxfold)
  ## draw heatmap A
  day_boxes <- data.frame(x=c(0.061, 0.6176),y=0.02)
  night_boxes <- data.frame(x=c(0.233, 0.791), y=0.02)
  heatA <- ggplot(figA.m, aes(variable, id_factor)) + 
           geom_tile(aes(fill = value)) +
           scale_fill_gradient2(name="Exp/Med", low = "blue", mid="grey20", high = "yellow", midpoint=1, space="Lab") +
           labs(title = mnames[1], x = "", y = "") + 
           scale_x_discrete(expand = c(0, 0)) +
           scale_y_discrete(expand = c(0, 0)) +
           theme(axis.ticks = element_blank(), 
                 plot.title = element_text(hjust = 0.5),
                 axis.line = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.text.y = element_blank() )
  ## draw heatmap B
  heatB <- ggplot(figB.m, aes(variable, id_factor)) + 
           geom_tile(aes(fill = value)) +
           scale_fill_gradient2(low = "blue", mid="grey20", high = "yellow", midpoint=1, space="Lab") + 
           labs(title = mnames[2], x = "", y = "") + 
           scale_x_discrete(expand = c(0, 0)) + 
           scale_y_discrete(expand = c(0, 0)) +
           theme(axis.ticks = element_blank(), 
                 plot.title = element_text(hjust = 0.5),
                 axis.line = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.text.y = element_blank(), 
                 legend.position="none")
  ## draw the two heatmap together and add light-dark box
  heat.out <- ggdraw() +
              draw_plot(heatA, 0.012, 0, 0.550, 0.998) +
              draw_plot(heatB, 0.568, 0, 0.405, 0.998) +
              geom_rect(data = day_boxes, aes(xmin = x, xmax = c(0.233, 0.791), ymin = y, ymax = y + 0.02), colour = "black", fill="white") +
              geom_rect(data = night_boxes, aes(xmin = x, xmax = c(0.399, 0.9595), ymin = y, ymax = y + 0.02), colour = "black", fill="black")
  ## output the heatmap
  save_plot(figfile, heat.out, ncol = 1, nrow = 1, base_height = 6.6, base_width = 6.1)
}
### read the file
exp_cirD <- read.csv("meta2d_seqA_liver_unigene.csv")
xr_cirD <- read.csv("meta2d_seqB_liver_unigene.csv")
tistag <- "liver"
minfold <- 0.65
maxfold <- 1.5
### draw the heatmap of circadian genes both datasets
bothD <- inner_join(exp_cirD, xr_cirD, by = "geneSym")
bothD <- data.frame (geneSym = bothD$geneSym)
fexp_cirD <- as.data.frame ( inner_join(bothD, exp_cirD, by = "geneSym") )  %>% 
  dplyr::select(geneSym, meta2d_pvalue, meta2d_phase, contains("GSM"))
fxr_cirD <- as.data.frame ( inner_join(bothD, xr_cirD, by = "geneSym") )  %>% 
  dplyr::select(geneSym, meta2d_pvalue, meta2d_phase, contains("ZT"))
### the required columns for input data frame: geneName, Pvalue, Phase, expression profile columns
figfile <- paste("compExp_", tistag, "_bothCir_heatmap.png", sep = "")
heatmapTCF(inputA=fexp_cirD, inputB=fxr_cirD, figfile = figfile, mnames = c("DataA", "DataB"), minfold=minfold, maxfold=maxfold)
