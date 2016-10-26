table <- read.table("gene_association_relevant_columns.csv", stringsAsFactors = F)

annTable <- NULL
for(i in unique(table[,1])) {
  rowNos <- grep(paste0("^", i, "$"), table[,1])
  temp <- c(i, paste(table[rowNos, 2], collapse = ","))
  annTable <- rbind(annTable, temp)
}
rownames(annTable) <- NULL
colnames(annTable) <- NULL

write.table(annTable, file = "ecoli_ann_table_for_topgo.txt", quote = F, sep = "\t", row.names = F, col.names = F)
