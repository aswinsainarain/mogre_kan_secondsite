# Load libraries
library(edgeR)

# Load geneReadCounts
geneReadCounts <- read.delim("../geneReadCounts.txt", stringsAsFactors = F)
rownames(geneReadCounts) <- geneReadCounts[,1]
geneReadCounts <- geneReadCounts[,-1]
colnames(geneReadCounts) <- c("cpxA_1", "cpxA_2", "cyaA_1", "cyaA_2", "fusA_1", "fusA_2",
                              "fusA-cpxA_1", "fusA-cpxA_2", "fusA-cyaA_1", "fusA-cyaA_2",
                              "fusA-rpoD_1", "fusA-rpoD_2", "fusA-topA_1", "fusA-topA_2",
                              "rpoD_1", "rpoD_2", "WT_1", "WT_2")

# make groups
group <- rep(c("cpxA", "cyaA", "fusA", "fusA-cpxA",
               "fusA-cyaA", "fusA-rpoD", "fusA-topA",
               "rpoD", "WT"), 
             each = 2)
dge <- DGEList(counts = geneReadCounts, group = group, genes = rownames(geneReadCounts))
dge <- calcNormFactors(dge)
write.table(dge$counts,
            file = "new_counts.txt",
            quote = F,
            sep = "\t",
            row.names = T,
            col.names = T)
## Get cpms
cpms <- cpm(dge)

## Get rpkms
# Getting gene lengths according to the row order in dge$counts
geneLengths <- read.table("../gene_lengths.txt", sep = "\t", stringsAsFactors = F)
geneLengths <- as.data.frame(geneLengths[,2], row.names  = geneLengths[,1])
colnames(geneLengths) <- "length"
geneLengths <- as.data.frame(geneLengths[rownames(dge$counts),]) # truncating the list according to the rows in the DGE object counts
rownames(geneLengths) <- rownames(dge$counts) # adopting the counts$counts rownames
colnames(geneLengths) <- "length"

rpkms <- rpkm(dge$counts, gene.length = as.vector(geneLengths[,1]))


# Filtering out genes that are not expressed in all samples
# Genes are kept if their cpms are at least one in any two samples, i.e the rowSum of cpms(dge)>1 (a T / F matrix) is >= 2
keep <- rowSums(cpms > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = F]
dge <- estimateDisp(dge)
#dge <- estimateCommonDisp(dge)
#dge <- estimateTagwiseDisp(dge)

# Perform the tests
strains <- levels(dge$samples$group)
strains <- strains[-length(strains)] ### removes the last entry i.e. WT

# the exactTestList will hold fold changes of all genes as determined by exactTest
exactTestList <- list()
for (i in 1:length(strains)) {
  exactTestList[[i]] <- exactTest(dge, pair = c("WT", strains[i]))
}
names(exactTestList) <- strains

exactTestTable <- NULL
for (i in 1:length(exactTestList)) {
  exactTestTable <- cbind(exactTestTable, exactTestList[[i]]$table$logFC)
}
colnames(exactTestTable) <- names(exactTestList)
rownames(exactTestTable) <- rownames(exactTestList[[1]]$table)
write.table(exactTestTable,
            file = "exactTestTable.txt",
            quote = F,
            sep = "\t",
            col.names = T,
            row.names = T)

# the decideTestsList will hold a vector of 0 and 1 / -1 indicating whether a gene is up / down or not mis regulated
decideTestsList <- lapply(exactTestList, decideTestsDGE, adjust.method = "BH", p.value = 0.001)
barplot(t(matrix(unlist(lapply(decideTestsList, summary)), nrow = 8, ncol = 3, byrow = T)[,-2]),
        names.arg = names(decideTestsList), 
        col = c("#2166AC","#B2182B"),
        border = NA)

summary(decideTestsDGE(exactTestList[[1]], adjust.method = "BH", p.value = 0.001))

## Whole list of all fold changes of all genes
wholeList <- list()
for (i in 1:length(exactTestList)) {
  wholeList[[i]] <- exactTestList[[i]]$table
}
names(wholeList) <- names(exactTestList)

## Toptags list of all fold changes of all genes with FDR
topTagsList <- list()
for (i in 1:length(exactTestList)) {
  topTagsList <- topTags(exactTestList[[i]], n = Inf, adjust.method = "BH") 
}

# sigList Only the genes that are de according to decideTest results; i.e. below p.value 0.001
# I'm not giving a fold change cutoff
sigList <- list()
for (i in 1:length(exactTestList)) {
  sigList[[i]] <- exactTestList[[i]]$table[(decideTestsList[[i]] == 1 | decideTestsList[[i]] == -1),]
}
names(sigList) <- names(exactTestList)

## Looking at fold changes of differntially expressed genes
par(mfrow = c(4,2)) 
for (i in 1:length(sigList) ) {
  hist(sigList[[i]]$logFC, breaks = 100,
       main = names(sigList)[i])
}
par(mfrow = c(4,2)) 
for (i in 1:length(sigList) ) {
  plot(density(sigList[[i]]$logFC),
       main = names(sigList)[i])
}
dev.off()

# What is the min fold change among up and down regulated genes
temp <- NULL
for (i in 1:length(sigList)) {
  temp <- c(temp, sigList[[i]]$logFC)
}
plot(density(temp))
min(temp[temp > 0]) # minimum up-regualted fold change
max(temp[temp < 0]) # minimum down-regulated fold change

# upsigList
upsigList <- list()
for (i in 1:length(sigList)) {
  upsigList[[i]] <- sigList[[i]][sigList[[i]]$logFC > 0, ]
}
names(upsigList) <- names(sigList)

# downsigList
downsigList <- list()
for (i in 1:length(sigList)) {
  downsigList[[i]] <- sigList[[i]][sigList[[i]]$logFC < 0, ]
}
names(downsigList) <- names(sigList)


# Getting percent of genes above log2 fold change of 1 or below -1 among sig list
sigListPerct <- NULL
for (i in 1:length(sigList)) {
  temp <- sigList[[i]]
  temp1 <- temp$logFC >= 1 | temp$logFC <= -1
  sigListPerct <- c(sigListPerct, (sum(temp1) / length(temp1) * 100)) 
}
names(sigListPerct) <- names(sigList)
