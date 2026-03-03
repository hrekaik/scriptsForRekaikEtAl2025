library(ggplot2)
library(ggrepel)

gitHubDir <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures/RNAseq/RNAseq_ActA/ActA_wtvsMutant"

deseq2.res <- read.delim(file.path(gitHubDir, "../../../outputs/RNAseq/RNAseq_ActA/ActA_wtvsMutant/DESeq2Analysis_96h.txt"))
hoxa.b.c.d.genes <- read.delim("../../../outputs/general/delBADC_genes.txt", header = FALSE)$V1
chr.y.genes <- deseq2.res$gene_short_name[grepl("chrY", deseq2.res$locus)]

deseq2.res$signif <- with(deseq2.res, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1.5)

deseq2.res$color <- "other"
deseq2.res$color[deseq2.res$gene_short_name %in% hoxa.b.c.d.genes] <- "HoxABCD"
deseq2.res$color[deseq2.res$gene_short_name %in% chr.y.genes] <- "chrY"

table(deseq2.res$signif)
# 855 genes

top50.genes <- deseq2.res[
  !is.na(deseq2.res$padj) &
    deseq2.res$padj < 0.05 &
    abs(deseq2.res$log2FoldChange) > 1.5,
]

top50.genes <- top50.genes[
  order(top50.genes$padj),
][1:min(50, nrow(top50.genes)), ]



deseq2.res$label <- ""
deseq2.res$label[deseq2.res$gene_short_name %in% top50.genes$gene_short_name] <-
  deseq2.res$gene_short_name[deseq2.res$gene_short_name %in% top50.genes$gene_short_name]


ggplot(deseq2.res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color)) +
  geom_text_repel(
    aes(label = label),
    box.padding = 0.5,
    max.overlaps = Inf
  ) +
  geom_point(
    data = subset(deseq2.res, label != ""),
    shape = 1,
    col = "red",
    size = 2
  ) +
  scale_color_manual(
    "Gene location",
    values = c("other" = "black", "HoxABCD" = "pink", "chrY" = "cyan")
  ) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  theme_classic()
ggsave("../../../outputs/RNAseq/RNAseq_ActA/ActA_wtvsMutant/Volcano_selectedgenes_96h.pdf", width = 7, height = 7)
