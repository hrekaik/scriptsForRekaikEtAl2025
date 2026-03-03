if (!"usefulLDfunctions" %in% installed.packages()) {
  devtools::install_github("lldelisle/usefulLDfunctions")
}
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("Signac")
safelyLoadAPackageInCRANorBioconductor("viridis")
safelyLoadAPackageInCRANorBioconductor("scCustomize")
safelyLoadAPackageInCRANorBioconductor("reshape2")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("ggseqlogo")
safelyLoadAPackageInCRANorBioconductor("ggrepel")

wd <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures"
plot.dir <- "outputs/Multiome/96h"

if (grepl("ldelisle", Sys.getenv("HOME"))) {
  wd <- "/scratch/ldelisle/rstudio_test/"
  plot.dir <- "scriptsForRekaikEtAl2025/outputs/Multiome/"
  dir.create(plot.dir, showWarnings = FALSE)
}
setwd(wd)

multiome.seurat <- readRDS("../toGEO/Multiome/96h/combined_Multiome_Velo_moreMotifs.RDS")
multiome.seurat$Genotype <- factor(x = multiome.seurat$Genotype, levels = c("wt", "BADC"))


# sample cells from wt and merge with BADC to do plots
set.seed(111)
BADC.cellnames <- WhichCells(multiome.seurat, expression = Genotype == "BADC")
wt.cellnames <- WhichCells(multiome.seurat, expression = Genotype == "wt")
wt.cellnames <- sample(wt.cellnames, size = length(BADC.cellnames), replace = F)
wt.BADC.cellnames <- c(wt.cellnames, BADC.cellnames)
multiome.seurat.subsampled <- subset(multiome.seurat, cells = wt.BADC.cellnames)

pdf(
  file = file.path(plot.dir, paste0("Multiome_bygeotype_wt.pdf")),
  width = 4, height = 4
)
DimPlot(multiome.seurat.subsampled,
  cols = c("#F5EDED", "#e69a8d"),
  label = T,
  pt.size = .1,
  alpha = 1,
  label.size = 4, group.by = "Genotype",
  order = c("BADC", "wt"),
  seed = 1
)
DimPlot(multiome.seurat.subsampled,
  cols = c("#F5EDED", "#5f4b8b"),
  label = T,
  pt.size = .1,
  alpha = 1,
  label.size = 4, group.by = "Genotype",
  order = c("wt", "BADC"),
  seed = 1
)
dev.off()


DefaultAssay(multiome.seurat.subsampled) <- "SCT"


my.markers <- c("T", "Tbx6", "Foxc1", "Pax3")
pdf(file = file.path(plot.dir, paste0("Multiome_geneexpression_Markers.pdf")),
    width = 16 , height = 4)
FeaturePlot_scCustom(multiome.seurat.subsampled , features = my.markers,
                     alpha_exp = 1,
                     pt.size = 1,
                     num_columns = length(my.markers))
dev.off()

my.markers <- c("Hoxa2", "Hoxa3", "Hoxb1")
pdf(file = file.path(plot.dir, paste0("Multiome_geneexpression_Hox.pdf")),
    width = 8 , height = 16)
FeaturePlot_scCustom(multiome.seurat.subsampled , features = my.markers,
                     alpha_exp = 1,
                     pt.size = 1,
                     split.by = "Genotype",
                     num_columns = 2)
dev.off()

my.markers <- c("Pbx1", "Meis1", "Pbx2", "Meis2")
pdf(file = file.path(plot.dir, paste0("Multiome_geneexpression_Pbx_Meis.pdf")),
    width = 8 , height = 16)
FeaturePlot_scCustom(multiome.seurat.subsampled , features = my.markers,
                     alpha_exp = 1,
                     pt.size = 1,
                     split.by = "Genotype",
                     num_columns = 2)
dev.off()



#### pseudo time clustering heatmap for ATAC peaks
caudal.tissue.seurat <- subset(multiome.seurat, cells = WhichCells(multiome.seurat,
                                                                   expression = seurat_clusters %in% 0:3))

# if (! file.exists("diffpeaks_allmainclusters.csv")) {
#   Variable_peaks_ATAC <- FindAllMarkers(
#     object = caudal.tissue.seurat,
#     assay = "ATAC",
#     min.pct = 0.1,
#     test.use = "wilcox",
#     only.pos = T
#   )
#   write.csv(Variable_peaks_ATAC, "diffpeaks_allmainclusters.csv")
# } else {
#   Variable_peaks_ATAC <- read.csv("diffpeaks_allmainclusters.csv")
# }

caudal.tissue.seurat$cluster.genotype <- paste0(caudal.tissue.seurat$seurat_clusters, "_", caudal.tissue.seurat$Genotype)
Idents(caudal.tissue.seurat) <- "cluster.genotype"
DimPlot(caudal.tissue.seurat, group.by = "cluster.genotype")
ggsave(file.path(plot.dir, paste0("Multiome_caudal_clustering.pdf")), width = 7, height = 7)

table(caudal.tissue.seurat$cluster.genotype)

if (!file.exists(file.path(plot.dir, "diffpeaks_allmainclustergenos.csv"))) {
  Variable_peaks_ATAC <- FindAllMarkers(
    object = caudal.tissue.seurat,
    assay = "ATAC",
    min.pct = 0.1,
    test.use = "wilcox",
    only.pos = T
  )
  write.csv(Variable_peaks_ATAC, file.path(plot.dir, "diffpeaks_allmainclustergenos.csv"))
} else {
  Variable_peaks_ATAC <- read.csv(file.path(plot.dir, "diffpeaks_allmainclustergenos.csv"))
}

# Select the top 500 peaks for the comparable clusters/genotype
# And remove duplicates

my.region.order <- c("1_wt", "1_BADC", "0_wt", "2_BADC")
topMarkers <- Variable_peaks_ATAC[Variable_peaks_ATAC$cluster %in% my.region.order, ]
topMarkers <- topMarkers[topMarkers$p_val_adj < 0.01, ]
toppeaks <- topMarkers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_log2FC)
toppeaks <- subset(toppeaks, !gene %in% toppeaks$gene[duplicated(toppeaks$gene)])
table(toppeaks$cluster)
# 0_wt   1_wt   2_wt   3_wt 2_BADC 0_BADC 1_BADC 3_BADC
#  500    476      0    500    500      0    476      0 
# 0_wt 1_BADC   1_wt 2_BADC   3_wt 
#  500    476    476    500    500 
toppeaks <- toppeaks[order(factor(toppeaks$cluster, levels = my.region.order)), ]

# Sort cells by pseudo time:

order.index <- order(caudal.tissue.seurat$velocity_pseudotime, decreasing = F)
length(order.index)
# 7252
cell.order <- colnames(caudal.tissue.seurat)[order.index]

####
# Define a function to smooth the openness along the pseudo time
rollMean <- function(prof, scale) {
  if (length(prof) <= scale) {
    return(NULL)
  }
  cprof <- cumsum(prof)
  if (scale %% 2 == 0) {
    d <- scale / 2
    return(
      c(
        cprof[d:(2 * d)] / d:(2 * d),
        (cprof[(2 * d + 1):length(cprof)] - cprof[1:(length(cprof) - 2 * d)]) / (2 * d),
        (cprof[length(cprof)] - cprof[(length(cprof) - 2 * d + 1):(length(cprof) - d)]) / (2 * d - 1):d
      )
    )
  } else {
    d <- (scale - 1) / 2
    return(
      c(
        cprof[d:(2 * d)] / d:(2 * d),
        cprof[scale] / scale,
        (cprof[(scale + 1):length(cprof)] - cprof[1:(length(cprof) - scale)]) / scale,
        (cprof[length(cprof)] - cprof[(length(cprof) - 2 * d):(length(cprof) - (d + 1))]) / (2 * d):(d + 1)
      )
    )
  }
}
####
DefaultAssay(caudal.tissue.seurat) <- "ATAC"
split.group <- "cluster.genotype"
my.cell.order <- c("1_wt", "0_wt", "1_BADC", "2_BADC")
caudal.tissue.seurat.subset <- subset(caudal.tissue.seurat, idents = my.cell.order)
caudal.tissue.seurat.subset$cluster.genotype <- factor(
  caudal.tissue.seurat.subset$cluster.genotype,
  levels = my.cell.order
)
ATAC.coverage.split.sorted <- lapply(
  unique(caudal.tissue.seurat.subset[[split.group]][, 1]),
  function(my.group.name) {
    FetchData(
      caudal.tissue.seurat.subset,
      cells = intersect(
        cell.order,
        colnames(caudal.tissue.seurat.subset)[caudal.tissue.seurat.subset[[split.group]][, 1] == my.group.name]
      ),
      vars = toppeaks$gene
    )
  }
)
ATAC.coverage.split.sorted.smooth <- lapply(
  ATAC.coverage.split.sorted,
  function(mat) {
    apply(
      X = mat,
      MARGIN = 2,
      FUN = rollMean,
      scale = 100
    )
  }
)
rm(ATAC.coverage.split.sorted)
ATAC.coverage.smooth.gg <- do.call(
  rbind,
  lapply(
    ATAC.coverage.split.sorted.smooth[!sapply(ATAC.coverage.split.sorted.smooth, is.null)],
    melt
  )
)
rm(ATAC.coverage.split.sorted.smooth)
# Add metadata and levels for plotting
ATAC.coverage.smooth.gg <- cbind(
  ATAC.coverage.smooth.gg,
  caudal.tissue.seurat.subset[[]][as.character(ATAC.coverage.smooth.gg$Var1), c(split.group, "velocity_pseudotime", "caudal_velocity_pseudotime")]
)
ATAC.coverage.smooth.gg$Var1 <- factor(ATAC.coverage.smooth.gg$Var1, levels = cell.order)
ATAC.coverage.smooth.gg$Var2 <- factor(ATAC.coverage.smooth.gg$Var2, levels = rev(toppeaks$gene))
ATAC.coverage.smooth.gg[, split.group] <- factor(
  ATAC.coverage.smooth.gg[, split.group],
  levels = levels(caudal.tissue.seurat.subset[[split.group]][, 1])
)

peaks.p <- ggplot(ATAC.coverage.smooth.gg, aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colours = (cividis(64))) +
  facet_grid(
    reformulate(split.group, "."),
    scales = "free", space = "free", drop = T
  )
peaks.p.rendered <- ggplot_gtable(ggplot_build(peaks.p))

pseudotime.p <- ggplot(ATAC.coverage.smooth.gg, aes(x = Var1, y = "Pseudotime", fill = velocity_pseudotime)) +
  geom_raster() +
  scale_fill_gradientn(colours = (viridis(64))) +
  facet_grid(
    reformulate(split.group, "."),
    scales = "free", space = "free", drop = T
  ) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
pseudotime.p.rendered <- ggplot_gtable(ggplot_build(pseudotime.p))
pseudotime.p.rendered$widths <- peaks.p.rendered$widths

DefaultAssay(caudal.tissue.seurat.subset) <- "RNA"
my.genes <- c("T", "Tbx6", "Foxc1", "Pax3")
RNA.expression.split.sorted <- lapply(
  unique(caudal.tissue.seurat.subset[[split.group]][, 1]),
  function(my.group.name) {
    FetchData(
      caudal.tissue.seurat.subset,
      cells = intersect(
        cell.order,
        colnames(caudal.tissue.seurat.subset)[caudal.tissue.seurat.subset[[split.group]][, 1] == my.group.name]
      ),
      vars = my.genes
    )
  }
)
RNA.expression.split.sorted.smooth <- lapply(
  RNA.expression.split.sorted,
  function(mat) {
    apply(
      X = mat,
      MARGIN = 2,
      FUN = rollMean,
      scale = 100
    )
  }
)
rm(RNA.expression.split.sorted)
RNA.expression.smooth.gg <- do.call(
  rbind,
  lapply(
    RNA.expression.split.sorted.smooth[!sapply(RNA.expression.split.sorted.smooth, is.null)],
    melt
  )
)
rm(RNA.expression.split.sorted.smooth)
# Add metadata and levels for plotting
RNA.expression.smooth.gg <- cbind(
  RNA.expression.smooth.gg,
  caudal.tissue.seurat.subset[[]][as.character(RNA.expression.smooth.gg$Var1), c(split.group, "velocity_pseudotime")]
)
RNA.expression.smooth.gg$Var1 <- factor(RNA.expression.smooth.gg$Var1, levels = cell.order)
RNA.expression.smooth.gg$Var2 <- factor(RNA.expression.smooth.gg$Var2, levels = rev(my.genes))
RNA.expression.smooth.gg[, split.group] <- factor(
  RNA.expression.smooth.gg[, split.group],
  levels = levels(caudal.tissue.seurat.subset[[split.group]][, 1])
)

rna.p <- ggplot(RNA.expression.smooth.gg, aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colours = (plasma(64))) +
  facet_grid(
    reformulate(split.group, "."),
    scales = "free", space = "free", drop = T
  ) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
rna.p.rendered <- ggplot_gtable(ggplot_build(rna.p))
rna.p.rendered$widths <- peaks.p.rendered$widths


pdf(
  file = file.path(plot.dir, paste0("doHeatmap_Multiome_peaks.pdf")),
  width = 8, height = 6
)
ggarrange(pseudotime.p.rendered, rna.p.rendered, peaks.p.rendered, ncol = 1, nrow = 3, heights = c(0.5, 1, 4))
dev.off()

# Finding overrepresented motifs with homer

toppeaks$seqnames <- sapply(strsplit(toppeaks$gene, "-"), "[", 1)
toppeaks$start <- sapply(strsplit(toppeaks$gene, "-"), "[", 2)
toppeaks$end <- sapply(strsplit(toppeaks$gene, "-"), "[", 3)

toppeaks.gr <- GenomicRanges::makeGRangesFromDataFrame(toppeaks, keep.extra.columns = TRUE)
toppeaks.gr.split <- split(toppeaks.gr, toppeaks.gr$cluster)
lapply(names(toppeaks.gr.split), function(geno.cluster) {
  rtracklayer::export.bed(GenomicRanges::sort(toppeaks.gr.split[[geno.cluster]]), file.path(plot.dir, paste0(geno.cluster, "_for_Homer.bed")))
})

# Plot enrichement of selected motifs
motif_list <- c("HOXA1(Homeobox)/mES-Hoxa1-ChIP-Seq(SRP084292)/Homer",
                "Meis1(Homeobox)/MastCells-Meis1-ChIP-Seq(GSE48085)/Homer",
                "JASPAR2022-MA1113.1")
for (n in seq_along(motif_list)) {
  pdf(file = file.path(plot.dir, paste0("Multiome_motifenrichement_", n, ".pdf")),
      width = 10 , height = 4.5)
  g <- FeaturePlot_scCustom(
    multiome.seurat.subsampled,
    colors_use = viridis(64, direction = -1),
    split.by = "Genotype",
    features = motif_list[n],
    alpha_exp = 1,
    pt.size = 1,
    num_columns = 2)
  print(g)
  dev.off()
}


# Quantitative Markers
caudal.tissue.seurat.subset$cell_type <- "Mesoderm"
caudal.tissue.seurat.subset$cell_type[caudal.tissue.seurat.subset$seurat_clusters == 1] <- "Progenitors"

# Specific genes
hoxa.b.c.d.genes <- read.delim(file.path(wd, "scriptsForRekaikEtAl2025/outputs/general/delBADC_genes.txt"), header = FALSE)$V1
chry.genes <- read.delim(file.path(wd, "scriptsForRekaikEtAl2025/outputs/general/chrYgenes.txt"), header = FALSE)$V1


DefaultAssay(caudal.tissue.seurat.subset) <- "RNA"

markers.RNA <- list()

# Progenitors
markers.fate.wtvsBADC <- FindMarkers(caudal.tissue.seurat.subset,
                                     ident.1 = "1_wt", ident.2 = "1_BADC",
                                     logfc.threshold = 0, min.pct = 0.01)
markers.fate.wtvsBADC$name <- rownames(markers.fate.wtvsBADC)
markers.fate.wtvsBADC$loci <- "other"
markers.fate.wtvsBADC$loci[rownames(markers.fate.wtvsBADC) %in% hoxa.b.c.d.genes] <- "HoxABCD"
markers.fate.wtvsBADC$loci[rownames(markers.fate.wtvsBADC) %in% chry.genes] <- "chrY"

write.table(markers.fate.wtvsBADC, file.path(plot.dir, "progenitors_wtVSBADC.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

significant.genes = subset(markers.fate.wtvsBADC, p_val_adj < 0.01 & (avg_log2FC > 1.5 | avg_log2FC < -1.5))
pdf(file = file.path(plot.dir, paste0("Volcano_progenitors.wtvsBADC.pdf")),
    width = 10 , height = 10)
g <- ggplot(data=markers.fate.wtvsBADC, aes(x=-avg_log2FC, y=-log10(p_val_adj))) +
  geom_point() + 
  theme_classic() +
  geom_text_repel(data=significant.genes,
                  aes(-avg_log2FC,-log10(p_val_adj),label=name)) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")
print(g)
dev.off()
pdf(file = file.path(plot.dir, paste0("Volcano_progenitors.wtvsBADC.pdf")),
    width = 10 , height = 10)
g <- ggplot(data=markers.fate.wtvsBADC, aes(x=-avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(aes(color = loci)) + 
  theme_classic() +
  geom_text_repel(data=significant.genes,
                  aes(-avg_log2FC,-log10(p_val_adj),label=name)) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red") +
  scale_color_manual("Gene location", values = c("other" = "black", "HoxABCD" = "pink", "chrY" = "cyan"))
print(g)
dev.off()

markers.RNA[["Progenitors"]] <- markers.fate.wtvsBADC
markers.RNA[["Progenitors"]]$group1 <- "Progenitors_wt"
markers.RNA[["Progenitors"]]$group2 <- "Progenitors_BADC"

# Mesoderm
markers.fate.wtvsBADC <- FindMarkers(caudal.tissue.seurat.subset,
                                     ident.1 = "0_wt", ident.2 = "2_BADC",
                                     logfc.threshold = 0, min.pct = 0.01)
markers.fate.wtvsBADC$name <- rownames(markers.fate.wtvsBADC)
markers.fate.wtvsBADC$loci <- "other"
markers.fate.wtvsBADC$loci[rownames(markers.fate.wtvsBADC) %in% hoxa.b.c.d.genes] <- "HoxABCD"
markers.fate.wtvsBADC$loci[rownames(markers.fate.wtvsBADC) %in% chry.genes] <- "chrY"

write.table(markers.fate.wtvsBADC, file.path(plot.dir, "mesoderm_wtVSBADC.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

significant.genes = subset(markers.fate.wtvsBADC, p_val_adj < 0.01 & (avg_log2FC > 1.5 | avg_log2FC < -1.5))
pdf(file = file.path(plot.dir, paste0("Volcano_mesoderm.wtvsBADC.pdf")),
    width = 10 , height = 10)
g <- ggplot(data=markers.fate.wtvsBADC, aes(x=-avg_log2FC, y=-log10(p_val_adj))) +
  geom_point() + 
  theme_classic() +
  geom_text_repel(data=significant.genes,
                  aes(-avg_log2FC,-log10(p_val_adj),label=name)) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")
print(g)
dev.off()
pdf(file = file.path(plot.dir, paste0("Volcano_mesoderm.wtvsBADC.pdf")),
    width = 10 , height = 10)
g <- ggplot(data=markers.fate.wtvsBADC, aes(x=-avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(aes(color = loci)) + 
  theme_classic() +
  geom_text_repel(data=significant.genes,
                  aes(-avg_log2FC,-log10(p_val_adj),label=name)) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red") +
  scale_color_manual("Gene location", values = c("other" = "black", "HoxABCD" = "pink", "chrY" = "cyan"))
print(g)
dev.off()

markers.RNA[["Mesoderm"]] <- markers.fate.wtvsBADC
markers.RNA[["Mesoderm"]]$group1 <- "Mesoderm_wt"
markers.RNA[["Mesoderm"]]$group2 <- "Mesoderm_BADC"

# Progenitors vs Mesoderm
Idents(caudal.tissue.seurat.subset) <- "cell_type"
markers.fate.wtvsBADC <- FindMarkers(caudal.tissue.seurat.subset,
                                     ident.1 = "Progenitors", ident.2 = "Mesoderm",
                                     logfc.threshold = 0, min.pct = 0.01)
Idents(caudal.tissue.seurat.subset) <- "cluster.genotype"
markers.fate.wtvsBADC$name <- rownames(markers.fate.wtvsBADC)
markers.fate.wtvsBADC$loci <- "other"
markers.fate.wtvsBADC$loci[rownames(markers.fate.wtvsBADC) %in% hoxa.b.c.d.genes] <- "HoxABCD"
markers.fate.wtvsBADC$loci[rownames(markers.fate.wtvsBADC) %in% chry.genes] <- "chrY"

write.table(markers.fate.wtvsBADC, file.path(plot.dir, "progenitorsVSmesodermal.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

significant.genes = subset(markers.fate.wtvsBADC, p_val_adj < 0.01 & (avg_log2FC > 1.5 | avg_log2FC < -1.5))
pdf(file = file.path(plot.dir, paste0("Volcano_progenitorsVSmesoderm.pdf")),
    width = 10 , height = 10)
g <- ggplot(data=markers.fate.wtvsBADC, aes(x=-avg_log2FC, y=-log10(p_val_adj))) +
  geom_point() + 
  theme_classic() +
  geom_text_repel(data=significant.genes,
                  aes(-avg_log2FC,-log10(p_val_adj),label=name)) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")
print(g)
dev.off()
pdf(file = file.path(plot.dir, paste0("Volcano_progenitorsVSmesoderm.pdf")),
    width = 10 , height = 10)
g <- ggplot(data=markers.fate.wtvsBADC, aes(x=-avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(aes(color = loci)) + 
  theme_classic() +
  geom_text_repel(data=significant.genes,
                  aes(-avg_log2FC,-log10(p_val_adj),label=name)) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red") +
  scale_color_manual("Gene location", values = c("other" = "black", "HoxABCD" = "pink", "chrY" = "cyan"))
print(g)
dev.off()

markers.RNA[["Both"]] <- markers.fate.wtvsBADC
markers.RNA[["Both"]]$group1 <- "Progenitors"
markers.RNA[["Both"]]$group2 <- "Mesoderm"



DefaultAssay(caudal.tissue.seurat.subset) <- "chromvar"

markers.chromvar <- list()

# Progenitors
markers.fate.wtvsBADC <- FindMarkers(caudal.tissue.seurat.subset,
                                     ident.1 = "1_wt", ident.2 = "1_BADC",
                                     mean.fxn = rowMeans,
                                     fc.name = "avg_diff")
markers.fate.wtvsBADC$name <- rownames(markers.fate.wtvsBADC)

write.table(markers.fate.wtvsBADC, file.path(plot.dir, "progenitors_chromvar_wtVSBADC.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

significant.genes = subset(markers.fate.wtvsBADC, p_val_adj < 0.01 & (avg_diff > 0.5 | avg_diff < -0.5))
pdf(file = file.path(plot.dir, paste0("Volcano_chromvar_progenitors.wtvsBADC.pdf")),
    width = 10 , height = 10)
g <- ggplot(data=markers.fate.wtvsBADC, aes(x=-avg_diff, y=-log10(p_val_adj))) +
  geom_point() + 
  theme_classic() +
  geom_text_repel(data=significant.genes,
                  aes(label=name)) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")
print(g)
dev.off()

markers.chromvar[["Progenitors"]] <- markers.fate.wtvsBADC
markers.chromvar[["Progenitors"]]$group1 <- "Progenitors_wt"
markers.chromvar[["Progenitors"]]$group2 <- "Progenitors_BADC"

# Mesodermal
markers.fate.wtvsBADC <- FindMarkers(caudal.tissue.seurat.subset,
                                     ident.1 = "0_wt", ident.2 = "2_BADC",
                                     mean.fxn = rowMeans,
                                     fc.name = "avg_diff")
markers.fate.wtvsBADC$name <- rownames(markers.fate.wtvsBADC)

write.table(markers.fate.wtvsBADC, file.path(plot.dir, "mesoderm_chromvar_wtVSBADC.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

significant.genes = subset(markers.fate.wtvsBADC, p_val_adj < 0.01 & (avg_diff > 0.5 | avg_diff < -0.5))
pdf(file = file.path(plot.dir, paste0("Volcano_chromvar_mesoderm.wtvsBADC.pdf")),
    width = 10 , height = 10)
g <- ggplot(data=markers.fate.wtvsBADC, aes(x=-avg_diff, y=-log10(p_val_adj))) +
  geom_point() + 
  theme_classic() +
  geom_text_repel(data=significant.genes,
                  aes(label=name)) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")
print(g)
dev.off()


markers.chromvar[["Mesoderm"]] <- markers.fate.wtvsBADC
markers.chromvar[["Mesoderm"]]$group1 <- "Mesoderm_wt"
markers.chromvar[["Mesoderm"]]$group2 <- "Mesoderm_BADC"

# Progenitors vs Mesoderm
Idents(caudal.tissue.seurat.subset) <- "cell_type"
markers.fate.wtvsBADC <- FindMarkers(caudal.tissue.seurat.subset,
                                     ident.1 = "Progenitors", ident.2 = "Mesoderm",
                                     mean.fxn = rowMeans,
                                     fc.name = "avg_diff")
Idents(caudal.tissue.seurat.subset) <- "cluster.genotype"
markers.fate.wtvsBADC$name <- rownames(markers.fate.wtvsBADC)

write.table(markers.fate.wtvsBADC, file.path(plot.dir, "chromvar_progenitorsVSmesodermal.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

significant.genes = subset(markers.fate.wtvsBADC, p_val_adj < 0.01 & (avg_diff > 0.5 | avg_diff < -0.5))
pdf(file = file.path(plot.dir, paste0("Volcano_chromvar_progenitorsVSmesodermal.pdf")),
    width = 10 , height = 10)
g <- ggplot(data=markers.fate.wtvsBADC, aes(x=-avg_diff, y=-log10(p_val_adj))) +
  geom_point() + 
  theme_classic() +
  geom_text_repel(data=significant.genes,
                  aes(label=name)) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")
print(g)
dev.off()


markers.chromvar[["Both"]] <- markers.fate.wtvsBADC
markers.chromvar[["Both"]]$group1 <- "Progenitors"
markers.chromvar[["Both"]]$group2 <- "Mesoderm"





my.genes <- c("Hoxa2", "Hoxa3", "Hoxb1",
              "Pbx1", "Pbx2", "Meis1", "Meis2",
              "Nkx1-2", "Cdx2", "Cdx1",
              "Wnt3a", "Fgf3", "Fgf4", "Etv4")
DotPlot(
  caudal.tissue.seurat.subset, features = my.genes,
  idents = c("1_wt", "1_BADC", "0_wt", "2_BADC"),
  assay = "RNA",
  scale = FALSE
)+
  theme(axis.text.x = element_text(angle = 90))

ggsave(file.path(plot.dir, "DotPlot_genes.pdf"), width = 7, height = 5)

markers.selected <- do.call(
  rbind,
  lapply(markers.RNA, subset, subset = name %in% my.genes)
)

write.table(markers.selected, file.path(plot.dir, "rna_selected.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)




tmp <- PseudobulkExpression(
  object = caudal.tissue.seurat.subset,
  group.by = "cluster.genotype",
  # features = my.genes, # I want to normalize to all genes so I must not set features here
  method = "aggregate",
  layer = "counts",
  return.seurat = TRUE
)
matrix.average.norm <- tmp[["RNA"]]$data[my.genes, ]
print(matrix.average.norm)
