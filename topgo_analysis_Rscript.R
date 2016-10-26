library(topGO)

######## GENE ONTOLOGY ANALYSIS USING TOPGO ############################
# generate ontology table
ontology <- readLines("go_analysis/cog_files/go-basic.obo")
ontoIds <- ontology[grep("^id: ", ontology)]
ontoNames <- ontology[grep("^name: ", ontology)]
ontoIds <- unlist(strsplit(ontoIds, split = "id: "))
ontoIds <- ontoIds[-seq(1, length(ontoIds), 2)]
ontoIds <- ontoIds[-((length(ontoIds)-4):length(ontoIds))]
ontoNames <- unlist(strsplit(ontoNames, split = "name: "))
ontoNames <- ontoNames[-seq(1, length(ontoNames), 2)]
ontoNames <- ontoNames[-((length(ontoNames)-4):length(ontoNames))]
ontology <- cbind(ontoIds, ontoNames)
rm(ontoIds, ontoNames)
rownames(ontology) <- ontology[,1]

# Using topgo # http://avrilomics.blogspot.in/2015/07/using-topgo-to-test-for-go-term.html

# importing gene annotations
geneID2GO <- readMappings(file = "go_analysis/cog_files/ecoli_ann_table_for_topgo.txt")
geneUniverse <- names(geneID2GO)# defining the gene universe

# # Testing
# genesOfInterest <- rownames(down_sig_WTfrt_vs_cpxA) # defining the genes of interest
# geneList <- factor(as.integer(geneUniverse %in% genesOfInterest)) # an vector that tells you which genes of interest are in the gene universe
# names(geneList) <- geneUniverse
# # constructing topGOdata object
# myGOdata <- new("topGOdata", 
#                 description="My project", 
#                 ontology="BP", 
#                 allGenes=geneList,  
#                 annot = annFUN.gene2GO, 
#                 gene2GO = geneID2GO)
# # looking at significant genes
# str(sigGenes(myGOdata))
# numSigGenes(myGOdata)
# # enrichment tests # p-values not corrected for multiple testing
# resultFisher <- runTest(myGOdata, 
#                         algorithm = "classic", 
#                         statistic = "fisher") 
# nodes2extract <- sum(score(resultFisher) < 0.01) # get number of significant nodes (GO terms)
# geneData(resultFisher) # summary of significant genes
# # all significant results
# allRes <- GenTable(myGOdata, 
#                    classic = resultFisher, 
#                    orderBy = "resultFisher", 
#                    ranksOf = "classic", 
#                    topNodes = nodes2extract
# )
# 
# # visualize position of statistically significant GO terms in the GO hierarchy
# showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
# # makes a pdf file tGO_classic_5_all.pdf
# printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

#### automated TOPGO analysis for all samples
geneOntUpList <- list()
geneOntDownList <- list()
for (i in 1:length(upsigList)) {
  tempUp <- rownames(upsigList[[i]])
  tempDown <- rownames(downsigList[[i]])
  
  # topgo analysis with up-regulated genes
  genesOfInterest <- tempUp # defining the genes of interest
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest)) # an vector that tells you which genes of interest are in the gene universe
  names(geneList) <- geneUniverse
  # constructing topGOdata object
  myGOdata <- new("topGOdata", 
                  description="My project", 
                  ontology="BP", 
                  allGenes=geneList,  
                  annot = annFUN.gene2GO, 
                  gene2GO = geneID2GO)
  # enrichment tests # p-values not corrected for multiple testing
  resultFisher <- runTest(myGOdata, 
                          algorithm = "classic", 
                          statistic = "fisher") 
  nodes2extract <- sum(score(resultFisher) < 0.01) # get number of significant nodes (GO terms)
  # all significant results
  allRes <- GenTable(myGOdata, 
                     classic = resultFisher, 
                     orderBy = "resultFisher", 
                     ranksOf = "classic", 
                     topNodes = nodes2extract
  )
  allRes[,"Term"] <- ontology[allRes[,1],2] # replaces truncated ontology terms with proper full terms
  # getting the genes
  temp <- genesInTerm(myGOdata)
  tempGenes <- NULL
  for (j in allRes[,1]) {
    tempGenes <- c(tempGenes, paste(temp[[j]][temp[[j]] %in% tempUp], collapse = ","))
  }
  allRes <- cbind(allRes, tempGenes)
  colnames(allRes) <- c(colnames(allRes)[-7], "genes")
  # writes the table to a file
  write.table(allRes,
              file = paste0("go_analysis/topgo_output/upregulated_genes/","topgo_up_", names(upsigList)[[i]], ".txt"),
              quote = F,
              sep = "\t",
              row.names = F,
              col.names = T)
  geneOntUpList[[i]] <- allRes
  # topgo analysis with down-regulated genes
  genesOfInterest <- tempDown # defining the genes of interest
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest)) # an vector that tells you which genes of interest are in the gene universe
  names(geneList) <- geneUniverse
  # constructing topGOdata object
  myGOdata <- new("topGOdata", 
                  description="My project", 
                  ontology="BP", 
                  allGenes=geneList,  
                  annot = annFUN.gene2GO, 
                  gene2GO = geneID2GO)
  # enrichment tests # p-values not corrected for multiple testing
  resultFisher <- runTest(myGOdata, 
                          algorithm = "classic", 
                          statistic = "fisher") 
  nodes2extract <- sum(score(resultFisher) < 0.01) # get number of significant nodes (GO terms)
  # all significant results
  allRes <- GenTable(myGOdata, 
                     classic = resultFisher, 
                     orderBy = "resultFisher", 
                     ranksOf = "classic", 
                     topNodes = nodes2extract
  )
  allRes[,"Term"] <- ontology[allRes[,1],2] # replaces truncated ontology terms with proper full terms
  # getting the genes
  temp <- genesInTerm(myGOdata)
  tempGenes <- NULL
  for (j in allRes[,1]) {
    tempGenes <- c(tempGenes, paste(temp[[j]][temp[[j]] %in% tempDown], collapse = ","))
  }
  allRes <- cbind(allRes, tempGenes)
  colnames(allRes) <- c(colnames(allRes)[-7], "genes")
  # write the table to a file
  write.table(allRes,
              file = paste0("go_analysis/topgo_output/downregulated_genes/","topgo_down_", names(downsigList)[[i]], ".txt"),
              quote = F,
              sep = "\t",
              row.names = F,
              col.names = T)
  geneOntDownList[[i]] <- allRes
}

names(geneOntUpList) <- names(geneOntDownList) <- names(upsigList)


### allTermsPvalsUpGenes table 
### contains Pvalues from topGO GO analysis for all up-regulated genes
allTerms <- NULL
for (i in 1:length(geneOntUpList)) {
  allTerms <- c(allTerms, geneOntUpList[[i]]$Term)
}
allTerms <- unique(allTerms)
allTerms <- sort(allTerms)

allTermsPvalsUpGenes <- NULL
for (i in 1:length(geneOntUpList)) {
  temp <- NULL
  for (j in 1:length(allTerms)) {
    if (allTerms[j] %in% geneOntUpList[[i]]$Term) {
      temp1 <- grep(paste0("^", allTerms[j], "$"), 
                    geneOntUpList[[i]]$Term)
      temp <- c(temp, geneOntUpList[[i]]$classic[temp1])
    } else {
      temp <- c(temp, "NA")
    }
  }
  allTermsPvalsUpGenes <- cbind(allTermsPvalsUpGenes, temp)
}
rownames(allTermsPvalsUpGenes) <- allTerms
colnames(allTermsPvalsUpGenes) <- names(geneOntUpList)
allTermsPvalsUpGenes <- apply(allTermsPvalsUpGenes, 1:2, as.numeric) # gives warnings - that's ok
allTermsPvalsUpGenes2 <- allTermsPvalsUpGenes
allTermsPvalsUpGenes2[is.na(allTermsPvalsUpGenes2)] <- 1

### allTermsPvalsDownGenes table 
### contains Pvalues from topGO GO analysis for all down-regulated genes
allTerms <- NULL
for (i in 1:length(geneOntDownList)) {
  allTerms <- c(allTerms, geneOntDownList[[i]]$Term)
}
allTerms <- unique(allTerms)
allTerms <- sort(allTerms)

### There are some problematic entries for grep among the down-regulated genes
### So I've generated independent if statements to take care of those entries
### phew!
allTermsPvalsDownGenes <- NULL
for (i in 1:length(geneOntDownList)) {
  temp <- NULL
  for (j in 1:length(allTerms)) {
    if (allTerms[j] == "[2Fe-2S] cluster assembly") {
      if (length(grep("^\\[2Fe-2S\\] cluster assembly$", geneOntDownList[[i]]$Term))!=0) {
        temp1 <- grep("^\\[2Fe-2S\\] cluster assembly$", geneOntDownList[[i]]$Term)
        if (length(temp1) > 1) { print (j)}
        temp <- c(temp, geneOntDownList[[i]]$classic[temp1])
      } else {
        temp <- c(temp, "NA")
      }
      
    } else if (allTerms[j] == "RNA (guanine-N7)-methylation") {
      if (length(grep("^RNA", allTerms[j], geneOntDownList[[i]]$Term)) != 0) {
        temp1 <- grep("^RNA", allTerms[j], geneOntDownList[[i]]$Term)
        if (length(temp1) > 1) { print (j)}
        temp <- c(temp, geneOntDownList[[i]]$classic[temp1])
      } else {
        temp <- c(temp, "NA")
      }
    } else if (allTerms[j] == "rRNA (guanine-N7)-methylation") {
      if (length(grep("^rRNA", allTerms[j], geneOntDownList[[i]]$Term)) != 0) {
        temp1 <- grep("^rRNA", allTerms[j], geneOntDownList[[i]]$Term)
        if (length(temp1) > 1) { print (j)}
        temp <- c(temp, geneOntDownList[[i]]$classic[temp1])
      } else {
        temp <- c(temp, "NA")
      }
    } else if (allTerms[j] %in% geneOntDownList[[i]]$Term) {
      temp1 <- grep(paste0("^", allTerms[j], "$"), 
                    geneOntDownList[[i]]$Term)
      if (length(temp1) > 1) { print (j)}
      temp <- c(temp, geneOntDownList[[i]]$classic[temp1])
    } else {
      temp <- c(temp, "NA")
    }
  }
  allTermsPvalsDownGenes <- cbind(allTermsPvalsDownGenes, temp)
}
rownames(allTermsPvalsDownGenes) <- allTerms
colnames(allTermsPvalsDownGenes) <- names(geneOntDownList)
allTermsPvalsDownGenes <- apply(allTermsPvalsDownGenes, 1:2, as.numeric) # gives warnings - that's ok
allTermsPvalsDownGenes2 <- allTermsPvalsDownGenes
allTermsPvalsDownGenes2[is.na(allTermsPvalsDownGenes2)] <- 1


#######################################################################

### Heatmaps of pvalues from topgo analysis
library(gplots)
library(fBasics)
## up-regulated genes
pdf("go_analysis/topgo_output/topgo_pvalues_heatmap_upregulated_genes.pdf", width = 10, height = 70)
heatmap.2(as.matrix(abs(log10(allTermsPvalsUpGenes2))),
          Rowv = T,
          Colv = T,
          col = rev(divPalette(100, name = "Spectral")),
          na.color = "grey",
          trace = "none",
          tracecol = NA,
          key.title = "Abs Log10 P-Value",
          main = "Up-regulated genes",
          cexCol = 1,
          cexRow = 0.95,
          margins = c(5,20),
          key.ylab = NA
          #key.ytickfun = NA
)
dev.off()
## down-regulated genes
pdf("go_analysis/topgo_output/topgo_pvalues_heatmap_downregulated_genes.pdf", width = 10, height = 80)
heatmap.2(as.matrix(abs(log10(allTermsPvalsDownGenes2))),
          Rowv = T,
          Colv = T,
          col = rev(divPalette(100, name = "Spectral")),
          na.color = "grey",
          trace = "none",
          tracecol = NA,
          key.title = "Abs Log10 P Value",
          main = "Down-regulated genes",
          cexCol = 1,
          cexRow = 0.95,
          margins = c(5,20),
          key.ylab = NA)
dev.off()