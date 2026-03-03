library(pheatmap)
directories <- c("RNAseq/RNAseq_EScells/", "RNAseq/RNAseq_Endo_diff/")
fpkms <- Reduce(merge, lapply(directories, function(dir){read.delim(file.path(dir, "outputs", "AllCufflinks_Simplified.txt.gz"))}))
samples.plans <- Reduce(rbind, lapply(directories, function(dir){read.delim(list.files(dir, "samplesPlan", full.names = TRUE))}))
rownames(samples.plans) <- samples.plans$sample
samples.plans$Replicate <- factor(samples.plans$Replicate)
samples.plans$Genotype <- samples.plans$Line
my.genes <- c("T", "Gsc", "Sox17", "Foxa2")
mat <- fpkms[match(my.genes, fpkms$gene_short_name), paste0("FPKM_", samples.plans$sample)]
colnames(mat) <- samples.plans$sample
rownames(mat) <- my.genes
clu <- hclust(dist(t(mat)))
dd <- as.dendrogram(clu)
clu2 <- reorder(dd, 1:nrow(samples.plans), agglo.FUN = mean)
pheatmap(
  log2(1 + mat),
  cluster_rows = FALSE,
  cluster_cols = as.hclust(clu2),
  show_colnames = FALSE,
  annotation_col = samples.plans[, c("Replicate", "Genotype", "Tissue")],
  annotation_colors = list("Genotype" = c("BADC" = "#e69a8d", "wt" = "#5f4b8b"),
                           "Replicate" = c("1" = "white", "2" = "grey70", "3" = "grey50"),
                           "Tissue" = c("Endo" = "orange", "ES_cells" = "grey30")),
  file = "outputs/RNAseq/RNAseq_both_combined/heatmap_diff.pdf",
  width = 5, height = 3
)
