####
# IMPORTANT FUNCTIONS!
# function to make a contingency table for overlapping gene sets (2) and perform the fischer test 
# Needs gene set A, B and a total number of objects
doFisherTest <- function(A, B, Total) {
  perOvlpAwithB <- (length(intersect(A, B)) / length(A)) * 100
  perOvlpBwithA <- (length(intersect(A, B)) / length(B)) * 100
  notAnotB <- Total - (length(A) + length(B) - length(intersect(A, B)))
  inAnotB <- length(A) - length(intersect(A, B))
  notAinB <- length(B) - length(intersect(A, B))
  inAinB <- length(intersect(A, B))
  overlapMatrix <- rbind(c(notAnotB, inAnotB),
                         c(notAinB, inAinB))
  rownames(overlapMatrix) <- c("notB", "inB")
  colnames(overlapMatrix) <- c("notA", "inA")
  return (list (percent_A_that_overlap_with_B = perOvlpAwithB,
                percent_B_that_overlap_with_A = perOvlpBwithA,
                overlapMatrix = overlapMatrix,
                fisher = fisher.test(overlapMatrix, alternative = "greater")))
}

# function to visualise all pairwise overlaps. 
# Give a list of all gene sets and the total no of genes and 
# whether to plot heat map with odds ratio (heatOdds = T) or with p-value (heatOdds = F)
# colPalette is the colours to be used for the heat map
doPairwiseFisher <- function (Sets, 
                              Total, 
                              heatOdds = T, 
                              colPalette = qualiPalette(100),
                              cellValueCol = "black") {
  library("gplots")
  library("fBasics")
  library("scales")  
  
  pValue <- NULL
  oddsRatio <- NULL
  
  for (i in 1:length(Sets)) {
    tempPValue <- NULL
    tempOddsRatio <- NULL
    for (j in 1:length(Sets)) {
      tempPValue <- c(tempPValue, doFisherTest(Sets[[i]], Sets[[j]], Total)$fisher$p.value)
      tempOddsRatio <- c(tempOddsRatio, doFisherTest(Sets[[i]], Sets[[j]], Total)$fisher$estimate)
    }
    pValue <- rbind(pValue, tempPValue)
    oddsRatio <- rbind(oddsRatio, tempOddsRatio)
  }
  
  pValue <- abs(log10(pValue))
  
  rownames(pValue) <- colnames(pValue) <- names(Sets)
  rownames(oddsRatio) <- colnames(oddsRatio) <- names(Sets)
  diag(pValue) <- NA
  diag(oddsRatio) <- NA
  
  if (heatOdds == T) {
    heatmap.2(oddsRatio, 
              symm = T,
              Rowv = F,
              Colv = F,
              col = colPalette, 
              na.color = "black", 
              trace = "none",
              density.info = "none",
              main = "Odds Ratio",
              cexRow = 1,
              cexCol = 1,
              cellnote = round(pValue, digits = 2), # Mutually exclusive to next line. Include if log10 p-value calculated.
              # cellnote = scientific(pValue, digits = 1), # Mutually exclusive to previous line. Include if log10 p-value is not calculated.
              notecol = cellValueCol,
              notecex = 1,
              mar = c(7, 12),
              key.title = "Odds Ratio",
              key.xlab = "")
    mtext("Cell Values: Abs Log10 p-Value", side = 3)
  } else {
    heatmap.2(pValue, 
              symm = T,
              Rowv = F,
              Colv = F,
              col = colPalette, 
              na.color = "black", 
              trace = "none",
              density.info = "none",
              main = "Abs log10 p-value",
              cexRow = 1,
              cexCol = 1,
              cellnote = round(oddsRatio, digits = 2), # Mutually exclusive to next line. Include if log10 p-value calculated.
              # cellnote = scientific(pValue, digits = 1), # Mutually exclusive to previous line. Include if log10 p-value is not calculated.
              notecol = cellValueCol,
              notecex = 1,
              mar = c(7, 12),
              key.title = "Abs Log10 p-value",
              key.xlab = "")
    mtext("Cell Values: Odds Ratio", side = 3)
  }
}

# function to visualise percentage of all pairwise overlaps. 
# Give a list of all sets
# colPalette is the colours to be used for the heat map - give only palette name in ""
# numCols is the number of columns in the plot as set by par(mfrow)
plotPercentOverlap <- function (Sets, 
                                colPalette = "qualiPalette",
                                cellValueCol = "black",
                                numCols = 2) {
  library("gplots")
  library("fBasics")
  
  colPalette = get(colPalette)(length(Sets))
  numRows = ceiling(length(Sets) / numCols)
  par(mfrow = c(numRows, numCols))
  par(mar = c(3,3,2,1))
  
  for (i in 1:length(Sets)) {
    A <- Sets[[i]]
    tempOvlp <- NULL
    tempOvlpText <- NULL
    for (j in 1:length(Sets)) {
      B <- Sets[[j]]
      tempOvlp <- c(tempOvlp, (length(intersect(A,B)) / length(A)) * 100)
      tempOvlpText <- c(tempOvlpText, length(intersect(A,B)))
    }
    barpos <- barplot(tempOvlp, 
                      names.arg = names(Sets), 
                      main = paste0(names(Sets[i]), ": ", length(A)), 
                      ylab = NULL, 
                      col = colPalette,
                      cex.names = 0.75,
                      #las = 3  # To rotate the x labels so that they are vertical instead of horizontal
    )
    #     if (i == 1) {  # In case a legend is required
    #       legend("top", 
    #              legend = names(Sets), 
    #              fill = colPalette, 
    #              bty = "n",
    #              cex = 0.85,
    #              #y.intersp = 0.5,
    #              ncol = 2)
    #     }
    text(x = barpos,
         y = 10,
         labels = tempOvlpText)
    mtext("% Overlap", side = 2, line = 2, cex = 0.65)
    
  }
}

# function to visualise percentage of overlap of the first list with remaining lists in the set 
# Give a list of all sets
# colPalette is the colours to be used for the heat map - give only palette name in ""
# numCols is the number of columns in the plot as set by par(mfrow)
plotPercentOverlap2 <- function (Sets, 
                                 colPalette = "qualiPalette",
                                 cellValueCol = "black",
                                 numCols = 2) {
  library("gplots")
  library("fBasics")
  par(mar = c(3,3,2,1))
  colPalette = get(colPalette)(length(Sets))
  A <- Sets[[1]]
  tempOvlp <- NULL
  tempOvlpText <- NULL
  for (j in 1:length(Sets)) {
    B <- Sets[[j]]
    tempOvlp <- c(tempOvlp, (length(intersect(A,B)) / length(A)) * 100)
    tempOvlpText <- c(tempOvlpText, length(intersect(A,B)))
  }
  barpos <- barplot(tempOvlp, 
                    names.arg = names(Sets), 
                    main = paste0(names(Sets[1]), ": ", length(A)), 
                    ylab = NULL, 
                    col = colPalette,
                    cex.names = 0.75,
                    #las = 3  # To rotate the x labels so that they are vertical instead of horizontal
  )
  #     if (i == 1) {  # In case a legend is required
  #       legend("top", 
  #              legend = names(Sets), 
  #              fill = colPalette, 
  #              bty = "n",
  #              cex = 0.85,
  #              #y.intersp = 0.5,
  #              ncol = 2)
  #     }
  text(x = barpos,
       y = 10,
       labels = tempOvlpText)
  mtext("% Overlap", side = 2, line = 2) 
}
####
### Function to plot venn diagrams given a named list x ###
plotMeAVenn <- function(x) {
  genes <- unique(as.vector(unlist(x)))
  VennMatrix <- NULL
  for (i in 1:length(x)) {
    VennMatrix <- cbind(VennMatrix, genes %in% x[[i]])
  }
  colnames(VennMatrix) <- names(x)
  library("venneuler")
  par(mar = c(1,1,1,1))
  plot(venneuler::venneuler(VennMatrix))
}
###
### Function to tell me the extent of all (binary) overlaps to fill in venn diagram
tellMeAllOverlaps <- function(x) {
  genes <- unique(as.vector(unlist(x)))
  overlapMatrix <- NULL
  for (i in 1:length(x)) {
    temp <- NULL
    for (j in 1:length(x)) {
      temp <- c(temp, length(intersect(x[[i]], x[[j]])))
    }
    overlapMatrix <- rbind(overlapMatrix, temp)
  }
  colnames(overlapMatrix) <- names(x)
  rownames(overlapMatrix) <- names(x)
  return(overlapMatrix)
}
##### List all functions #####
lsf.str()
##############################
