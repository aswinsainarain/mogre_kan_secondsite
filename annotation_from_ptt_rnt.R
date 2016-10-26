### Gene annotations wherever possible ##################################################################
ptt <- read.delim("NC_000913.ptt", skip = 2)  # skip first two lines
rnt <- read.delim("NC_000913.rnt", skip = 2)  # skip first two lines
annotation <- NULL
for (i in 1:nrow(ptt)) {
  temp <- unlist(strsplit(as.character(ptt[i,1]), split = "[..]")) # 1 is start 3 is stop
  temp1 <- NULL
  temp1 <- c(as.character(ptt[i,5]), temp[1], temp[3])
  annotation <- rbind(annotation, temp1)
}
for (i in 1:nrow(rnt)) {
  temp <- unlist(strsplit(as.character(rnt[i,1]), split = "[..]")) # 1 is start 3 is stop
  temp1 <- NULL
  temp1 <- c(as.character(rnt[i,5]), temp[1], temp[3])
  annotation <- rbind(annotation, temp1)
}
annotation <- cbind(as.data.frame(annotation[,1]), apply(annotation[,2:3], 2, as.numeric ))
rownames(annotation) <- NULL
annotation <- annotation[order(annotation[,2]),]
annotation <- annotation[-(grep("^ins", annotation[,1])),] # remove all ins
write.table(annotation, file = "mg1655.annotation", quote = F, sep = "\t", row.names = F, col.names = c("Gene", "Start", "Stop"))
##########################################################################################################