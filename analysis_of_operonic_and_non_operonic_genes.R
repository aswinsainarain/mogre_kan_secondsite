# annotation file
# ann <- read.delim("mg1655_strand_info_added.annotation", stringsAsFactors = F)
# igRegions <- NULL
# for (i in 2:nrow(ann)) {
#   a <- paste0(ann[(i-1), 1], "-", ann[i, 1])
#   b <- ann[(i-1), 3] + 1
#   c <- ann[i, 2] - 1
#   igRegions <- rbind(igRegions, c(a, b, c))
# }
# rownames(igRegions) <- igRegions[,1]
# igRegions <- igRegions[,-1]
# igRegions <- apply(igRegions, 1:2, as.numeric)
# igRegions <- igRegions[!((igRegions[,2] - igRegions[,1]) <= 0),]
# igRegions <- cbind(igRegions, (igRegions[,2] - igRegions[,1]))
# colnames(igRegions) <- c("Start", "End", "Length")
# igRegions <- as.data.frame(igRegions)
# 
# perNtCts <- read.table("../wtfrt_1.readsPerNt")
# rownames(perNtCts) <- as.character(perNtCts[,1])

# This is point-less. It is very slow.
# igRegionsCounts <- NULL
# for (i in 1:nrow(igRegions)) {
#   j <- (igRegions[i, 2]:igRegions[i, 3]) 
#   j <- as.character(j)
#   igRegionsCounts <- c(igRegionsCounts, sum(perNtCts[j, 2]))
# }


### Operon analysis ###
geneReadCounts <- read.delim("../geneReadCounts.txt", stringsAsFactors = F)
rownames(geneReadCounts) <- geneReadCounts[,1]
geneReadCounts <- geneReadCounts[,-1]
colnames(geneReadCounts) <- c("cpxA_1", "cpxA_2", "cyaA_1", "cyaA_2",
                              "fusA_1", "fusA_2", "fusA-cpxA_1", "fusA-cpxA_2",
                              "fusA-cyaA_1", "fusA-cyaA_2", "fusA-rpoD_1", "fusA-rpoD_2",
                              "fusA-topA_1", "fusA-topA_2", "rpoD_1", "rpoD_2", "WT_1", "WT_2")

ptt <- read.delim("NC_000913.ptt", skip = 2, stringsAsFactors = F)

# operon list obtained from DOOR
# http://csbl.bmb.uga.edu/DOOR/displayNCoperon.php?id=1944&page=1&nc=NC_000913#tabs-1
door <- read.delim("door_e_coli_mg1655_operon_list.opr", stringsAsFactors = F)
rownames(ptt) <- ptt$Synonym
door <- cbind(door, ptt[door$Synonym, "Gene"])
colnames(door)[10] <- "Gene"
door[,10] <- as.character(door[,10])
door$OperonID <- as.character(door$OperonID)

operons <- unique(door$OperonID)
# how many genes per operon so that I take only those that have >= 2 genes i.e. real operons
genesPerOperon <- NULL
for (i in 1:length(operons)) {
  temp <- grep(paste0("^",operons[i], "$"), door$OperonID)
  genesPerOperon <- rbind(genesPerOperon, c(operons[i], length(temp)))
}
genesPerOperon <- as.data.frame(genesPerOperon, stringsAsFactors = F)
genesPerOperon[,2] <- as.numeric(genesPerOperon[,2])
boxplot(genesPerOperon[,2])

realOperons <- genesPerOperon[genesPerOperon[,2] > 1, ] # these guys even define operons with a single gene
boxplot(realOperons[,2])

nonOperonicGenes <- NULL
# first get the genes that are in operons
for (i in 1:nrow(realOperons)) {
  temp <- grep(paste0("^", realOperons[i,1], "$"), door$OperonID)
  nonOperonicGenes <- c(nonOperonicGenes, door$Gene[temp])
}
# then take the genes that are not in operons
nonOperonicGenes <- ptt$Gene[!(ptt$Gene %in% nonOperonicGenes)]

# Plot of SDs of read counts in operonic genes and non-operonic genes
# I'm calculating SDs of each operon separately. Hence there will be a distribution of SDs.
# On the other hand there is just a single list of non operonic genes. Thus there will be a single median in the boxplot.
# Also although it is incorrect, sd has been calculated for just two genes if present in an operon
# Function to find sample SDs i.e. divided by n and not n-1
sampleSD <- function (temp) {
  toReturn <- sqrt(  sum((temp - mean(temp, na.rm = T))^2, na.rm = T) / length(temp))
  return(toReturn)
}

## Plotting sample SDs

pdf("operonic_vs_non_opernoic.pdf", width = 5, height = 20)
par(mfrow = c(9,2), mar = c(4,4,2,1))
for (j in 1:ncol(geneReadCounts)) {
  
  operonReadCounts <- list()
  for (i in 1:nrow(realOperons)) {
    temp <- grep(paste0("^", realOperons[i,1], "$"), door$OperonID)
    temp <- door$Gene[temp]
    operonReadCounts[[i]] <- geneReadCounts[temp, j]
  }
  
  operonicSds <- unlist(lapply(operonReadCounts, sampleSD))
  boxplot(operonicSds,
          sampleSD(geneReadCounts[nonOperonicGenes, j]),
          #ylim = c(0, 3000),
          names = c("Operonic", "Non Operonic"),
          ylab = "SD of read counts",
          main = colnames(geneReadCounts)[j]
  )
  
}
dev.off()

pdf("operonic_vs_non_opernoic_with_ylim.pdf", width = 5, height = 20)
par(mfrow = c(9,2), mar = c(4,4,2,1))
for (j in 1:ncol(geneReadCounts)) {
  
  operonReadCounts <- list()
  for (i in 1:nrow(realOperons)) {
    temp <- grep(paste0("^", realOperons[i,1], "$"), door$OperonID)
    temp <- door$Gene[temp]
    operonReadCounts[[i]] <- geneReadCounts[temp, j]
  }
  
  operonicSds <- unlist(lapply(operonReadCounts, sampleSD))
  boxplot(operonicSds,
          sampleSD(geneReadCounts[nonOperonicGenes, j]),
          ylim = c(0, 3000),
          names = c("Operonic", "Non Operonic"),
          ylab = "SD of read counts",
          main = colnames(geneReadCounts)[j]
  )
  
}
dev.off()


#### Aswin's suggestion ####
# Each operonic SD can be divided by SD of all genes
# This should be below 1
# You can even take a random set of genes from the wildtype and calculate their SDs and plot their distribution
# This time I'm dividing the read counts by length of each gene as well since the gene lengths in an operon could vary
# I'm not dividing by library size since that will not change 

geneLengths <- NULL
for(i in 1:nrow(ptt)) {
  temp <- ptt[i, 1]
  temp <- as.numeric(unlist(strsplit(temp, split = "[.][.]")))
  temp <- temp[2] - temp[1]
  geneLengths <- c(geneLengths, temp)
}
names(geneLengths) <- ptt$Gene
geneLengths <- geneLengths[rownames(geneReadCounts)]

operonicSdsTable <- NULL
for (j in 1:ncol(geneReadCounts)) {
  
  operonReadCounts <- list()
  for (i in 1:nrow(realOperons)) {
    temp <- grep(paste0("^", realOperons[i,1], "$"), door$OperonID)
    temp <- door$Gene[temp]
    operonReadCounts[[i]] <- geneReadCounts[temp, j] / geneLengths[temp]
  }
  
  operonicSds <- unlist(lapply(operonReadCounts, sampleSD))
  operonicSdsTable <- cbind(operonicSdsTable, operonicSds)

}
colnames(operonicSdsTable) <- colnames(geneReadCounts)
rownames(operonicSdsTable) <- rownames(realOperons)

# getting normalized gene read counts to get genomic average SD
normGRC <- geneReadCounts / geneLengths
sampleAvgSd <- apply(normGRC, 2, sampleSD)

# operonic SDs divided by average sample sd # all normalised read counts
# the transposing is essential otherwise you will get incorrect results
opByAvgSds <- t(t(operonicSdsTable) / sampleAvgSd)
colnames(opByAvgSds) <- colnames(geneReadCounts)
rownames(opByAvgSds) <- rownames(realOperons)
boxplot(opByAvgSds, las = 2 
        #ylim = c(0, 0.1)
)

# Now I have to randomly sample genes equal to the number of genes in each operon and get their sd
# Will have to do this for each sample
randomTable <- NULL
for (i in 1:ncol(normGRC)) {
  temp <- NULL
  for (j in 1:nrow(realOperons)) {
    temp2 <- sample(1:nrow(normGRC), realOperons[j,2]) # sample row numbers 
    temp2 <- normGRC[temp2, i]
    temp <- c(temp, (sampleSD(temp2) / sampleSD(normGRC[,i])))
  }
  randomTable <- cbind(randomTable, temp)
}
colnames(randomTable) <- colnames(normGRC)
rownames(randomTable) <- rownames(realOperons)
boxplot(randomTable, las = 2)

# Are the means of the real and random different
pVals <- NULL
for (i in 1:ncol(opByAvgSds)) {
  pVals <- c(pVals, wilcox.test(opByAvgSds[,i], randomTable[,i])$p.value)
}
names(pVals) <- colnames(opByAvgSds)
max(pVals)

# boxplots of real and random
pdf("real_vs_random_sd_per_op_div_by_total_sd.pdf", width = 10, height = 6)
par(mfrow = c(2,1), mar = c(4,4,2,1))
boxplot(opByAvgSds, las = 2, main = "Real", cex.axis = 0.6, ylab = "Per op SD / Total SD")
boxplot(randomTable, las = 2, main = "Random", cex.axis = 0.6, ylab = "Per op SD / Total SD")
dev.off()
pdf("real_vs_random_sd_per_op_div_by_total_sd_with_ylim.pdf", width = 10, height = 6)
par(mfrow = c(2,1), mar = c(5.2,4,2,1))
boxplot(opByAvgSds, las = 2, main = "Real", cex.axis = 0.8, ylim = c(0, 0.2), notch = T, ylab = "Per op SD / Total SD")
boxplot(randomTable, las = 2, main = "Random", cex.axis = 0.8, ylim = c(0, 0.2), notch = T, ylab = "Per op SD / Total SD")
dev.off()