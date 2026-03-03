library(ggplot2)
library(ggrepel)
deseq2.res <- read.delim("outputs/RNAseq/RNAseq_Endo_diff/DESeq2AnalysisForFactorLine.txt")
hoxa.b.c.d.genes <- read.delim("outputs/general/delBADC_genes.txt", header = FALSE)$V1
chr.y.genes <- deseq2.res$gene_short_name[grepl("chrY", deseq2.res$locus)]

deseq2.res$signif <- with(deseq2.res, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1.5)

deseq2.res$color <- "other"
deseq2.res$color[deseq2.res$gene_short_name %in% hoxa.b.c.d.genes] <- "HoxABCD"
deseq2.res$color[deseq2.res$gene_short_name %in% chr.y.genes] <- "chrY"

table(deseq2.res$signif)
# 855 genes

selected.genes <- read.delim("RNAseq/RNAseq_Endo_diff/geneslist.txt")
deseq2.res$label <- ""
deseq2.res$label[deseq2.res$gene_short_name %in% selected.genes$gene_short_name] <-
  deseq2.res$gene_short_name[deseq2.res$gene_short_name %in% selected.genes$gene_short_name]

ggplot(deseq2.res, aes(x = log2FoldChange, y = -log(padj))) +
  geom_point(aes(color = color)) +
  geom_text_repel(aes(label = label), box.padding = 0.5, max.overlaps = Inf) +
  geom_point(data = subset(deseq2.res, label != ""), shape = 1, col = "red") +
  scale_color_manual("Gene location", values = c("other" = "black", "HoxABCD" = "pink", "chrY" = "cyan")) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  theme_minimal()
ggsave("outputs/RNAseq/RNAseq_Endo_diff/Volcano_selectedgenes.pdf", width = 7, height = 7)
