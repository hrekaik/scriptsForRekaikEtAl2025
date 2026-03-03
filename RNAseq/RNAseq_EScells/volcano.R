library(ggplot2)
library(ggrepel)
deseq2.res <- read.delim("outputs/RNAseq/RNAseq_EScells/DESeq2AnalysisForFactorLine.txt")
hoxa.b.c.d.genes <- read.delim("outputs/general/delBADC_genes.txt", header = FALSE)$V1
chr.y.genes <- deseq2.res$gene_short_name[grepl("chrY", deseq2.res$locus)]

deseq2.res$signif <- with(deseq2.res, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1.5)

deseq2.res$color <- "other"
deseq2.res$color[deseq2.res$gene_short_name %in% hoxa.b.c.d.genes] <- "HoxABCD"
deseq2.res$color[deseq2.res$gene_short_name %in% chr.y.genes] <- "chrY"

table(deseq2.res$signif)
# 49 genes

ggplot(deseq2.res, aes(x = log2FoldChange, y = -log(padj))) +
  geom_point(aes(color = color)) +
  geom_text_repel(data = subset(deseq2.res, signif), aes(label = gene_short_name), size = 2) +
  scale_color_manual("Gene location", values = c("other" = "black", "HoxABCD" = "pink", "chrY" = "cyan")) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme_minimal()
ggsave("outputs/RNAseq/RNAseq_EScells/Volcano.pdf", width = 7, height = 7)
