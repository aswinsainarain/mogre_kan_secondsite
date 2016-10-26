# Load the edgeR workspace before using this script

# Network of transcription factors table from regulonDB
tfNetwork <- read.table("tf_target_analysis/tf_targets_lists/network_tf_gene.txt", 
                        sep = "\t", comment.char = "#", stringsAsFactors = F)
sigmaNetwork <- read.table("tf_target_analysis/tf_targets_lists/network_sigma_gene.txt", 
                           sep = "\t", comment.char = "#", stringsAsFactors = F)


##### Making a list that has all targets of all TFs ################
allGenesEcoli <- rownames(geneLengths)
allNetworksEcoli <- rbind(tfNetwork, sigmaNetwork)
allNetworksEcoli <- allNetworksEcoli[allNetworksEcoli[,2] %in% allGenesEcoli, ]
allRegulatorsEcoli <- unique(allNetworksEcoli[,1])

# A list of vectors. Each vector contains named + or - entries.
# Each vector denotes a TF. Each element in the vector denotes target of that TF.
# + or - indicates whether it is positively regulated or negatively regulated

tfTargetInteractionList <- list()
for (i in 1:length(unique(allNetworksEcoli[,1]))) {
  temp0 <- paste0("^", unique(allNetworksEcoli[,1])[i], "$")
  temp <- allNetworksEcoli[grep(temp0, allNetworksEcoli[,1]),]
  ##### geting + - or ? for all gene targets of the (i)th transcription factor
  temp2 <- NULL
  for (j in temp[duplicated(temp[,2]), 2]) {
    temp1 <- temp[grep(j, temp[,2]),3]
    if ("+" %in% temp1 & "-" %in% temp1 | "+-" %in% temp1 | "-+" %in% temp1) {
      temp2 <- rbind(temp2, c(j, "+-"))
    }
    if (("+" %in% temp1 & "-" %in% temp1 | "+-" %in% temp1 | "-+" %in% temp1) & "?" %in% temp1) {
      temp2 <- rbind(temp2, c(j, "+-?"))
    }
    if (("+" %in% temp1 & "?" %in% temp1) & !("-" %in% temp1)) {
      temp2 <- rbind(temp2, c(j, "+?"))
    }
    if (("-" %in% temp1 & "?" %in% temp1) & !("+" %in% temp1)) {
      temp2 <- rbind(temp2, c(j, "-?"))
    }
  }
  
  #for (j in temp2[duplicated(temp2[,1]), 1]) {
  #  print(temp2[grep(j, temp2[,1]),])
  #}
  # I've checked that duplicated elements in temp2 have the same second column and hence can be removed
  temp2 <- temp2[!(duplicated(temp2[,1])),]
  temp3 <- as.data.frame(temp[(!(temp[,2] %in% temp[duplicated(temp[,2]), 2])), 2:3])
  colnames(temp3) <- c("V1", "V2")
  temp2 <- rbind(temp3, temp2)
  rownames(temp2) <- temp2[,1]
  # ok so temp2 has + - or ? of all gene targets of the (i)th transcription factor
  temp3 <- temp2[,2]
  names(temp3) <- rownames(temp2)
  tfTargetInteractionList[[i]] <- temp3
}
names(tfTargetInteractionList) <- unique(allNetworksEcoli[,1])

temp <- NULL
for(i in tfTargetInteractionList) {
  temp <- c(temp, names(i))
}
temp <- unique(temp)
tfTargetInteractionTable <- matrix(data = NA, 
                                   nrow = length(temp), 
                                   ncol = length(tfTargetInteractionList))
tfTargetInteractionTable <- as.data.frame(tfTargetInteractionTable)
rownames(tfTargetInteractionTable) <- temp
colnames(tfTargetInteractionTable) <- names(tfTargetInteractionList)
for (i in 1:length(tfTargetInteractionList)) {
  temp1 <- temp %in% names(tfTargetInteractionList[[i]])
  tfTargetInteractionTable[temp1,i] <- as.vector(tfTargetInteractionList[[i]])
}
write.table(tfTargetInteractionTable, 
            file = "tf_target_analysis/tfTargetInteractionTable.txt",
            quote = F,
            sep = "\t",
            row.names = T,
            col.names = T)

#######################################################################################################################################

## Making a table that has the number of targets of all TFs among the up-regulated genes of all mutants
targetsInUp <- NULL
for(i in 1:length(upsigList)) {
  temp <- NULL
  for (j in 1:length(tfTargetInteractionList)) {
    temp1 <- names(tfTargetInteractionList[[j]]) %in% rownames(upsigList[[i]])
    temp1 <- sum(as.numeric(temp1))
    temp <- c(temp, temp1)
  }
  targetsInUp <- cbind(targetsInUp, temp)
}
rownames(targetsInUp) <- names(tfTargetInteractionList)
colnames(targetsInUp) <- names(upsigList)
targetsInUp2 <- targetsInUp[rowSums(targetsInUp) > 10,]

## Making a table that has the number of targets of all TFs among the down-regulated genes of all mutants
targetsInDown <- NULL
for(i in 1:length(downsigList)) {
  temp <- NULL
  for (j in 1:length(tfTargetInteractionList)) {
    temp1 <- names(tfTargetInteractionList[[j]]) %in% rownames(downsigList[[i]])
    temp1 <- sum(as.numeric(temp1))
    temp <- c(temp, temp1)
  }
  targetsInDown <- cbind(targetsInDown, temp)
}
rownames(targetsInDown) <- names(tfTargetInteractionList)
colnames(targetsInDown) <- names(downsigList)
targetsInDown2 <- targetsInDown[rowSums(targetsInDown) > 10,]

### Getting P values for targets in up and down
### Using Fisher Test
### Total number of genes ? = 4042 courtsey nrow(dge$counts)
# Need to run the functions_Rscript.R before running this part! 

# Generating P vals for the targetsInUp table 
targetsInUpPVals <- NULL
for (i in 1:ncol(targetsInUp)) {
  tempGeneSet1 <- rownames(upsigList[[colnames(targetsInUp)[i]]])
  tempColumn <- NULL
  for (j in 1:nrow(targetsInUp)) {
    tempGeneSet2 <- names(tfTargetInteractionList[[rownames(targetsInUp)[j]]])
    tempColumn <- c(tempColumn,
                    doFisherTest(A = tempGeneSet1, 
                                 B = tempGeneSet2, 
                                 Total = 4042)$fisher$p.value)
  }
  targetsInUpPVals <- cbind(targetsInUpPVals, tempColumn)
  rm(tempColumn, tempGeneSet1, tempGeneSet2)
}
colnames(targetsInUpPVals) <- colnames(targetsInUp)
rownames(targetsInUpPVals) <- rownames(targetsInUp)

# Generating P vals for the targetsInDown table 
targetsInDownPVals <- NULL
for (i in 1:ncol(targetsInDown)) {
  tempGeneSet1 <- rownames(downsigList[[colnames(targetsInDown)[i]]])
  tempColumn <- NULL
  for (j in 1:nrow(targetsInDown)) {
    tempGeneSet2 <- names(tfTargetInteractionList[[rownames(targetsInDown)[j]]])
    tempColumn <- c(tempColumn,
                    doFisherTest(A = tempGeneSet1, 
                                 B = tempGeneSet2, 
                                 Total = 4042)$fisher$p.value)
  }
  targetsInDownPVals <- cbind(targetsInDownPVals, tempColumn)
  rm(tempColumn, tempGeneSet1, tempGeneSet2)
}
colnames(targetsInDownPVals) <- colnames(targetsInDown)
rownames(targetsInDownPVals) <- rownames(targetsInDown)
#####






### counting + / - / other targets among the up and down genes in the mutants
# cyaA mutant
temp <- rownames(upsigList[["cyaA"]])
temp <- tfTargetInteractionTable[temp, "CRP"]
temp <- as.factor(temp)
temp <- summary(temp)

temp1 <- rownames(downsigList[["cyaA"]])
temp1 <- tfTargetInteractionTable[temp1, "CRP"]
temp1 <- as.factor(temp1)
temp1 <- summary(temp1)

temp2 <- rbind(temp1[1:3], temp[1:3])
rownames(temp2) <- c("Down", "Up")

write.table(temp2,
            "tf_target_analysis/cyaA_mutant_crp_targets_table.txt",
            quote = F,
            sep = "\t")

# fusA2-cyaA mutant
temp <- rownames(upsigList[["fusA2-cyaA"]])
temp <- tfTargetInteractionTable[temp, "CRP"]
temp <- as.factor(temp)
temp <- summary(temp)

temp1 <- rownames(downsigList[["fusA2-cyaA"]])
temp1 <- tfTargetInteractionTable[temp1, "CRP"]
temp1 <- as.factor(temp1)
temp1 <- summary(temp1)

temp2 <- rbind(temp1[1:2], temp[1:2])
rownames(temp2) <- c("Down", "Up")

write.table(temp2,
            "tf_target_analysis/fusA_cyaA_mutant_crp_targets_table.txt",
            quote = F,
            sep = "\t")


### counting + / - / other targets among the up and down genes in the cpxA mutants
sum(!is.na(tfTargetInteractionTable[,"CpxR"])) # Total number of targets of CpxR
# cpxA mutant
temp <- rownames(upsigList[["cpxA"]])
temp <- tfTargetInteractionTable[temp, "CpxR"]
temp <- as.factor(temp)
temp <- summary(temp)

temp1 <- rownames(downsigList[["cpxA"]])
temp1 <- tfTargetInteractionTable[temp1, "CpxR"]
temp1 <- as.factor(temp1)
temp1 <- summary(temp1)

temp2 <- rbind(temp1[1:2], temp[1:2])
rownames(temp2) <- c("Down", "Up")

write.table(temp2,
            "tf_target_analysis/cpxA_mutant_CpxR_targets_table.txt",
            quote = F,
            sep = "\t")

# fusA2-cpxA mutant
temp <- rownames(upsigList[["fusA2-cpxA"]])
temp <- tfTargetInteractionTable[temp, "CpxR"]
temp <- as.factor(temp)
temp <- summary(temp)

temp1 <- rownames(downsigList[["fusA2-cpxA"]])
temp1 <- tfTargetInteractionTable[temp1, "CpxR"]
temp1 <- as.factor(temp1)
temp1 <- summary(temp1)

temp2 <- rbind(temp1[1:2], c(temp[1], 0))
rownames(temp2) <- c("Down", "Up")

write.table(temp2,
            "tf_target_analysis/fusA_cpxA_mutant_CpxR_targets_table.txt",
            quote = F,
            sep = "\t")

## Getting list of common genes among the cpxA and cpxA.fusA mutants
temp <- unique(c(rownames(sigList[["cpxA"]]), rownames(sigList[["fusA2-cpxA"]])))
temp <- temp[temp %in% rownames(sigList[["cpxA"]]) & temp %in% rownames(sigList[["fusA2-cpxA"]])]
plot(sigList[["cpxA"]][temp, 1],
     sigList[["fusA2-cpxA"]][temp, 1],
     xlab = "cpxA",
     ylab = "fusA2-cpxA",
     xlim = c(-5,4),
     ylim = c(-5,4))
abline(a = 0, b = 1)

temp1 <- sigList[["cpxA"]][temp, 1] / sigList[["fusA2-cpxA"]][temp, 1]
names(temp1) <- temp
barplot(temp1, ylab = "cpxA / fusA2-cpxA", las = 2)
abline(h = 1)
sigList[["cpxA"]][c("cpxA", "cpxP", "cpxR"),]
sigList[["fusA2-cpxA"]][c("cpxA", "cpxP", "cpxR"),]

temp2 <- tfTargetInteractionTable[temp, "CpxR"]
tfTargetInteractionTable[c("cpxA", "cpxP", "cpxR"), "CpxR"] # all are positive targets


## Fold changes of siderophore metabolism genes across mutants
entGenes <- names(c(tfTargetInteractionList[["Sigma19"]],
                    tfTargetInteractionList[["Fur"]]))
entGenesFCs <- NULL
for (i in 1:length(wholeList)) {
  entGenesFCs <- cbind(entGenesFCs, wholeList[[i]][entGenes, "logFC"])
}
colnames(entGenesFCs) <- names(wholeList)
rownames(entGenesFCs) <- entGenes

library(gplots)
library(fBasics)
pdf("tf_target_analysis/fold_changes_of_targets_of_fur_and_sigma19.pdf",
    width = 5, height = 20)
heatmap.2(entGenesFCs,
          Rowv = T,
          Colv = T,
          col = rev(divPalette(100, name = "RdBu")),
          trace = "none",
          tracecol = NA,
          key.title = "log2 FC",
          margins = c(7,5)
          )
dev.off()