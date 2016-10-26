fileList <- list.files(pattern = ".perGeneFRCounts$")
png(filename = "fwd_vs_rvrse_plot.png", width = 1000, height = 1000)
par(mfrow = c(6,6))
#allCorr <- NULL
for (i in fileList) {
  data <- read.table(file = i, sep = "\t", header = T, stringsAsFactors = F)
  row.names(data) <- data[,1]
  data <- data[,-1]
  corValue <- cor(data[,1], data[,2])
  #allCorr <- c(allCorr, corValue)
  plot(data[,1], data[,2], xlab = "# forward", ylab = '# reverse', xlim = c(0, 1000), ylim = c(0,1000), main = i)
  text(800, 800, labels = paste("Pcorr:", round(corValue, digits = 3)))
}
dev.off()

pdf(file = "fwd_vs_rvrse_plot.pdf", width = 15, height = 15)
par(mfrow = c(6,6))
for (i in fileList) {
  data <- read.table(file = i, sep = "\t", header = T, stringsAsFactors = F)
  row.names(data) <- data[,1]
  data <- data[,-1]
  corValue <- cor(data[,1], data[,2])
  plot(data[,1], data[,2], xlab = "# forward", ylab = '# reverse', xlim = c(0, 1000), ylim = c(0,1000), main = i)
  text(800, 800, labels = paste("Pcorr:", round(corValue, digits = 3)))
}
dev.off()
