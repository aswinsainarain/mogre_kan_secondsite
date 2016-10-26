### don't save this workspace

### normalising the data myself using mode
### and then finding fold changes
### just an additional check
getMode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
apply(dge$counts, 2, getMode)
normCounts <- dge$counts / apply(dge$counts, 2, getMode)

meanNormCounts <- NULL
for (i in seq(1, ncol(normCounts), 2)) {
  meanNormCounts <- cbind(meanNormCounts, rowMeans(normCounts[,(i:(i+1))]))
}
colnames(meanNormCounts) <- c(strains, "WTfrt")
foldChanges <-  meanNormCounts / meanNormCounts[,"WTfrt"]
logFC <- apply(foldChanges, 1:2, log2)
logFC[is.na(logFC)] <- NA
logFC[is.infinite(logFC)] <- NA
pdf("edgeR_vs_manual_log2FC_scatter_plots.pdf", width = 6, height = 10)
par(mfrow = c(4,2), mar = c(4,4,2,1))
for (i in 1:ncol(exactTestTable)) {
  smoothScatter(exactTestTable[,i], logFC[,i],
                xlab = "edgeR log2FC", ylab = "manual log2FC",
                main = strains[i])
  legend("topleft",
         legend = paste("Cor:", round(cor(exactTestTable[,i], logFC[,i], use = "complete.obs"), 3)),
         bty = "n")
}
dev.off()