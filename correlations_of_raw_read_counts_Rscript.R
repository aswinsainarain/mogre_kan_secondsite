# needs EdgeR workspace
library(gplots)
library(fBasics)

temp <- cor(geneReadCounts)

pdf("correlations_of_raw_read_counts.pdf", width = 8, height = 7)
heatmap.2(as.matrix(temp),
          Rowv = F,
          Colv = F,
          breaks = seq(0, 1, length.out = 101),
          col = divPalette(100, name = "Spectral"),
          trace = "none",
          tracecol = NA,
          key.title = NA,
          key.xlab = "Correlation",
          margins = c(7,7))
dev.off()
