
options(stringsAsFactors = FALSE)
rm(list = ls())

library(ggplot2)

gitHubDir <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures/RNAseq/RNAseq_ActA/ActA_wtvsMutant"


# Expression table (genes x samples)
tableWithNormalizedExpression <- file.path(gitHubDir, "../../../outputs/RNAseq/RNAseq_ActA/ActA_wtvsMutant/AllCufflinks_Simplified.txt")

# Samples plan
samplesPlan <- file.path(gitHubDir,"samplesPlan_ActA_wtvsBADC.txt")

# Use FPKM columns?
useFPKM <- TRUE

# Number of most variable genes to keep
restrictToNMoreVariantGenes <- 1000

# Output file
outputFile <- file.path(gitHubDir,"../../../outputs/RNAseq/RNAseq_ActA/ActA_wtvsMutant/PCA_PC1_PC2.pdf")


expr <- read.delim(tableWithNormalizedExpression, check.names = FALSE)
samplesPlanDF <- read.delim(samplesPlan)

if (!"sample" %in% colnames(samplesPlanDF)) {
  stop("samplesPlan must contain a column named 'sample'")
}

if (useFPKM) {
  colnames(expr) <- gsub("^FPKM_", "", colnames(expr))
}

samplesToPlot <- intersect(colnames(expr), samplesPlanDF$sample)

if (length(samplesToPlot) < 2) {
  stop("Less than 2 samples found in common between expression and samplesPlan")
}


exprData <- expr[, samplesToPlot]

# Remove genes with zero expression everywhere
exprData <- exprData[rowSums(exprData) != 0, ]

# Log-transform
exprData <- log2(exprData + 1)


if (!is.null(restrictToNMoreVariantGenes)) {
  vars <- apply(exprData, 1, var)
  exprData <- exprData[
    order(vars, decreasing = TRUE)[
      seq_len(min(nrow(exprData), restrictToNMoreVariantGenes))
    ],
  ]
}


pca <- prcomp(t(exprData), center = TRUE, scale. = FALSE)

# Percent variance
varExplained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)


pcaDF <- data.frame(
  sample = rownames(pca$x),
  pca$x
)

pcaDF <- merge(pcaDF, samplesPlanDF, by = "sample")


#PCA PLOT (PC1 vs PC2)
p <- ggplot(
  pcaDF,
  aes(
    x = PC1,
    y = PC2,
    color = Line,
    shape = Tissue
  )
) +
  geom_point(size = 6) +
  theme_classic(base_size = 16) +
  scale_color_manual(values = c( "#e69a8d","#5f4b8b")) +
  xlab(paste0("PC1 (", varExplained[1], "%)")) +
  ylab(paste0("PC2 (", varExplained[2], "%)")) +
  scale_x_continuous(limits = c(-55, 55)) +
  scale_y_continuous(limits = c(-30, 40)) +
  labs(title = "PCA of RNA-seq samples")

ggsave(outputFile, p, width = 7, height = 6)

cat("PCA plot saved to:", outputFile, "\n")
