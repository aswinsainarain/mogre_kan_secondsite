library(gplots)
library(fBasics)

########### Reading Peter et. al. dataset ############################################
# peters supercoiling sensitive genes only
peters <-  readLines("peter_et_al_supplementary/peters_gene_exp_ratios_of_supercoiling_sensitive_genes.csv")
peters <- peters[-1]
temp <- NULL
for (i in peters) {
  temp <- rbind(temp, strsplit(i, split = "\t")[[1]])
}
peters <- temp[-1,-1]
tempNames <- unlist(strsplit(temp[,1], split = " "))[-1]
for (i in 1:length(tempNames)) {
  if (tempNames[i] %in% ptt[, "Synonym"]) {
    tempNames[i] <- ptt[grep(tempNames[i], ptt[, "Synonym"]), "Gene"]
  }
}
rownames(peters) <- tempNames
colnames(peters) <- temp[1,][-1]
rm(temp)
rm(tempNames)
peters <- as.data.frame(peters) 
peters <- cbind(peters[,1:2], apply(peters[,3:ncol(peters)], 2, as.numeric))
# "correl_sigma" and "p-value"
# 106 were induced (+ve correlation) by plasmid relaxation and 
# 200 were repressed (-ve correlation) by plasmid relaxation.
plot(density(peters[,"correl_sigma"]))
ssg <- rownames(peters)
ssg <- ssg[ssg %in% rownames(ptt)]

# peters all genes
petersAll <-  readLines("peter_et_al_supplementary/peters_exp_ratios_of_all_genes.csv")
petersAll <- petersAll[-1]
temp <- NULL
for (i in petersAll) {
  temp <- rbind(temp, strsplit(i, split = "\t")[[1]])
}
petersAll <- temp[-1,-1]
tempNames <- unlist(strsplit(temp[,1], split = " "))[-1]
for (i in 1:length(tempNames)) {
  if (tempNames[i] %in% ptt[, "Synonym"]) {
    tempNames[i] <- ptt[grep(tempNames[i], ptt[, "Synonym"]), "Gene"]
  }
}
tempNames <- tempNames[-c(338,492,941,3494)]
petersAll <- petersAll[-c(338,492,941,3494),]
rownames(petersAll) <- tempNames
colnames(petersAll) <- temp[1,][-1]
rm(temp)
rm(tempNames)
petersAll <- as.data.frame(petersAll) 
petersAll <- cbind(petersAll[,1:2], apply(petersAll[,3:ncol(petersAll)], 2, as.numeric))


#################### Reading datasets done ###########################################


### Checking expression of topA, topB, parC, parE, gyrA, gyrB, in Novo and Norflox data
pdf("peter_et_al_supplementary/exp_ratios_of_supercoiling_genes_in_petersAll.pdf", 
    width = 8, height = 4.5)
supGenes <- petersAll[c("topA", "topB", "parC", "parE", "gyrA", "gyrB"), ]
heatmap.2(as.matrix(supGenes[,-(1:5)]),
          Rowv = F,
          Colv = F,
          col = rev(divPalette(100, name = "RdBu")),
          trace = "none",
          tracecol = NA,
          key.title = "Exp ratios",
          margins = c(6,5),
          colsep = c(5, 10, 14, 18, 27, 32)
          )
dev.off()

#
plot(wholeList[["fusA.topA"]][,1], petersAll[rownames(wholeList[["fusA.topA"]]), "nov200_1"])
cor(wholeList[["fusA.topA"]][,1], petersAll[rownames(wholeList[["fusA.topA"]]), "nov200_1"],
    use = "complete", method = "spearman")

plot(wholeList[["fusA.topA"]][,1], petersAll[rownames(wholeList[["fusA.topA"]]), "nor20wt"])
cor(wholeList[["fusA.topA"]][,1], petersAll[rownames(wholeList[["fusA.topA"]]), "nor20wt"],
    use = "complete", method = "spearman")

temp <- unique(rownames(wholeList[["fusA.topA"]]), rownames(petersAll))
#

# Overlap with known supercoiling sensitive genes
pdf("peter_et_al_supplementary/ssg_overlap_with_de_genes_topA_mutant.pdf", width = 5, height = 3)
temp <- list (SSG = ssg,
              fusA.topA = rownames(sigList[["fusA.topA"]]))
temp <- list (SSG = ssg,
              fusA.topA_up = rownames(upsigList[["fusA.topA"]]),
              fusA.topA_down = rownames(downsigList[["fusA.topA"]]))
plotMeAVenn(temp)
dev.off()
tellMeAllOverlaps(temp)
doPairwiseFisher(Sets = temp,
                 Total = nrow(dge$counts),
                 heatOdds = F)
###


# The supercoiling correlation is based on their entire datset. 
# Can I look at correlations of my dataset with any of theirs perhaps?
pdf("peters_all_conditions_corr.pdf", width = 7, height = 5.5)
#png("peters_all_conditions_corr.png", width = 700, height = 550)
tempNames <- names(petersAll)
genes <- unique(c(rownames(petersAll), rownames(wholeList[["fusA.topA"]])))
genes <- genes[genes %in% rownames(petersAll) & genes %in% rownames(wholeList[["fusA.topA"]])]
PwithPetersAllCorrMatrix <- cbind(wholeList[["fusA.topA"]][genes, "logFC"],
                                  petersAll[genes,6:40])
PwithPetersAllCorrMatrix <- cor(PwithPetersAllCorrMatrix, method = "spearman")
rownames(PwithPetersAllCorrMatrix) <- colnames(PwithPetersAllCorrMatrix) <- c("fusA.topA", colnames(petersAll[,6:40]))
pdf("peter_et_al_supplementary/topA_mutant_correlation_with_peter_dataset.pdf",
    width = 8, height = 8)
heatmap.2(apply(as.matrix(PwithPetersAllCorrMatrix), c(1,2), as.numeric),
          Rowv = F,
          Colv = F,
          col = rev(divPalette(100, name = "RdBu")),
          colsep = c(1,6,11,15,19,28,33),
          rowsep = c(1,6,11,15,19,28,33),
          trace = "none",
          tracecol = NA,
          cexRow = 0.9,
          cexCol = 0.9,
          margins = c(6,6),
          key.title = "Correlation"
)
dev.off()



# "correl_supercoiling" and "p-value"
plot(density(petersAll[,"correl_supercoiling"]))
# Now comparison with my dataset
# Scatter plots of supercoiling correlation with log2FCs
# Correlation bet
par(mfrow = c(2,3))
par(mar = c(4,4,1,1))
# all up-regulated genes in WTfrt_0K_vs_P_0K
genes <- upExpressed$up_sig_WTfrt_0K_vs_P_0K[upExpressed$up_sig_WTfrt_0K_vs_P_0K %in% rownames(petersAll)]
plot(petersAll[genes, "correl_supercoiling"], WTfrt_0K_vs_P_0K$table[genes, 1],
     xlab = "Peter et. al. supercoiling corr",
     ylab = "P_0K up Log2FC")
temp <- round(cor(petersAll[genes, "correl_supercoiling"], WTfrt_0K_vs_P_0K$table[genes, 1]), 3)
legend("topleft", legend = paste("Corr", temp), bty = "n", text.col = "blue")
# all down-regulated genes in WTfrt_0K_vs_P_0K
genes <- downExpressed$down_sig_WTfrt_0K_vs_P_0K[downExpressed$down_sig_WTfrt_0K_vs_P_0K %in% rownames(petersAll)]
plot(petersAll[genes, "correl_supercoiling"], WTfrt_0K_vs_P_0K$table[genes, 1],
     xlab = "Peter et. al. supercoiling corr",
     ylab = "P_0K down Log2FC")
temp <- round(cor(petersAll[genes, "correl_supercoiling"], WTfrt_0K_vs_P_0K$table[genes, 1]), 3)
legend("bottomleft", legend = paste("Corr", temp), bty = "n", text.col = "blue")

# all genes
genes <- rownames(WTfrt_0K_vs_P_0K)[rownames(WTfrt_0K_vs_P_0K) %in% rownames(petersAll)]
plot(petersAll[genes, "correl_supercoiling"], WTfrt_0K_vs_P_0K$table[genes, 1],
     xlab = "Peter et. al. supercoiling corr",
     ylab = "P_0K Log2FC")
temp <- round(cor(petersAll[genes, "correl_supercoiling"], WTfrt_0K_vs_P_0K$table[genes, 1]), 3)
legend("topleft", legend = paste("Corr", temp), bty = "n", text.col = "blue")

# all up-regulated genes in WTfrt_0K_vs_P_1K
genes <- upExpressed$up_sig_WTfrt_0K_vs_P_1K[upExpressed$up_sig_WTfrt_0K_vs_P_1K %in% rownames(petersAll)]
plot(petersAll[genes, "correl_supercoiling"], WTfrt_0K_vs_P_1K$table[genes, 1],
     xlab = "Peter et. al. supercoiling corr",
     ylab = "P_1K up Log2FC")
temp <- round(cor(petersAll[genes, "correl_supercoiling"], WTfrt_0K_vs_P_1K$table[genes, 1]), 3)
legend("topleft", legend = paste("Corr", temp), bty = "n", text.col = "blue")

# all down-regulated genes in WTfrt_0K_vs_P_1K
genes <- downExpressed$down_sig_WTfrt_0K_vs_P_1K[downExpressed$down_sig_WTfrt_0K_vs_P_1K %in% rownames(petersAll)]
plot(petersAll[genes, "correl_supercoiling"], WTfrt_0K_vs_P_1K$table[genes, 1],
     xlab = "Peter et. al. supercoiling corr",
     ylab = "P_1K down Log2FC")
temp <- round(cor(petersAll[genes, "correl_supercoiling"], WTfrt_0K_vs_P_1K$table[genes, 1]), 3)
legend("bottomleft", legend = paste("Corr", temp), bty = "n", text.col = "blue")

# all genes
genes <- rownames(WTfrt_0K_vs_P_1K)[rownames(WTfrt_0K_vs_P_1K) %in% rownames(petersAll)]
plot(petersAll[genes, "correl_supercoiling"], WTfrt_0K_vs_P_1K$table[genes, 1],
     xlab = "Peter et. al. supercoiling corr",
     ylab = "P_1K Log2FC")
temp <- round(cor(petersAll[genes, "correl_supercoiling"], WTfrt_0K_vs_P_1K$table[genes, 1]), 3)
legend("topleft", legend = paste("Corr", temp), bty = "n", text.col = "blue")

dev.off()

# The supercoiling correlation is based on their entire datset. 
# Can I look at correlations of my dataset with any of theirs perhaps?
pdf("peters_all_conditions_corr.pdf", width = 7, height = 5.5)
#png("peters_all_conditions_corr.png", width = 700, height = 550)
tempNames <- names(petersAll)
genes <- rownames(WTfrt_0K_vs_P_0K)[rownames(WTfrt_0K_vs_P_0K) %in% rownames(petersAll)]
PwithPetersAllCorrMatrix <- cbind(WTfrt_0K_vs_P_0K$table[genes,1],
                                  WTfrt_0K_vs_P_1K$table[genes,1],
                                  petersAll[genes,6:40])
PwithPetersAllCorrMatrix <- cor(PwithPetersAllCorrMatrix)
rownames(PwithPetersAllCorrMatrix) <- colnames(PwithPetersAllCorrMatrix) <- c("FusA.TopA 0K", "FusA.TopA.1K",
                                                                              colnames(petersAll[,6:40]))
heatmap.2(apply(as.matrix(PwithPetersAllCorrMatrix), c(1,2), as.numeric),
          Rowv = F,
          Colv = F,
          col = rev(divPalette(100, name = "RdBu")),
          trace = "none",
          cexRow = 0.75,
          cexCol = 0.75,
          margins = c(6,6),
          main = "Correlation"
)
rm(genes)
rm(tempNames)
dev.off()

# Overlap with known supercoiling sensitive genes
pdf("ssg_overlap_with_de_genes.pdf", width = 5, height = 3)
temp <- list (SSG = ssg,
              DE_genes = unique(c(diffExpressed$sig_WTfrt_0K_vs_P_0K,
                                  diffExpressed$sig_WTfrt_0K_vs_P_1K)))
plotMeAVenn(temp)
dev.off()
tellMeAllOverlaps(temp)
doFisherTest(temp[[1]], temp[[2]], Total = Total)

## check fold expression of ssg vs non-ssg in topA mutants
pdf("kdensity_ssg_non-ssg_in_topA_strains.pdf", width = 7, height = 3)
de_ssg <- intersect(temp[[1]], temp[[2]])
de_non_ssg <- setdiff(temp[[2]], temp[[1]])

colNames <- c("firebrick1", "dodgerblue1")
par(mfrow = c(1,2))
par(mar = c(4,4,2,2))
plot(density(WTfrt_0K_vs_P_0K$table[de_ssg, 1]), col = colNames[1], main = "P_0K", lwd = 2,
     xlab = expression("Log"[2]*" fold change"), ylim =c(0, 0.35))
points(density(WTfrt_0K_vs_P_0K$table[de_non_ssg, 1]), type = "l", col = colNames[2], lwd = 2)
legend("topleft", col = colNames, lty = 1, lwd = 2, legend = c("ssg", "non ssg"), bty = "n")
plot(density(WTfrt_0K_vs_P_1K$table[de_ssg, 1]), col = colNames[1], main = "P_1K", lwd = 2,
     xlab = expression("Log"[2]*" fold change"), ylim = c(0, 0.35))
points(density(WTfrt_0K_vs_P_1K$table[de_non_ssg, 1]), type = "l", col = colNames[2], lwd = 2)
legend("topleft", col = colNames, lty = 1, lwd = 2, legend = c("ssg", "non ssg"), bty = "n")
dev.off()
