library(usefulLDfunctions)

safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("scCustomize")
safelyLoadAPackageInCRANorBioconductor("tidyr")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("viridis")
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("pheatmap")
safelyLoadAPackageInCRANorBioconductor("ggrepel")


wd <- "/Papers/DelHox/toGEO/scRNAseq"
plot.dir <- "../../Code_figures/outputs/scRNAseq/"


setwd(wd)

nameRDS <- "combined_AlltimewtvsBADC_96f83823.analyzed.scVelo.RDS"
combined.seurat <- readRDS(nameRDS)


pdf(file = file.path(plot.dir, "bytime.pdf"),
    width = 5 , height = 4)
DimPlot(combined.seurat,
        cols = c("#9CA986","#C9DABF","#5F6F65","#423E28"), 
        label = T,
        pt.size = 1,
        alpha = 1,
        label.size = 4, group.by = "Time",
        shuffle = T,
        seed = 1)
dev.off()

pdf(
  file = file.path(plot.dir, paste0("bygenotype_BADC.pdf")),
  width = 6, height = 4
)
DimPlot(combined.seurat,
        cols = c("#F5EDED", "#e69a8d"),
        label = T,
        pt.size = .1,
        alpha = 1,
        label.size = 4, group.by = "Genotype",
        order = c("BADC", "wt"),
        seed = 1
)
dev.off()

table(combined.seurat$Genotype,combined.seurat$Time)

pdf(
  file = file.path(plot.dir, paste0("bygenotype_wt.pdf")),
  width = 6, height = 4
)
DimPlot(combined.seurat,
        cols = c("#F5EDED", "#5f4b8b"),
        label = T,
        pt.size = .1,
        alpha = 1,
        label.size = 4, group.by = "Genotype",
        order = c("wt", "BADC"),
        seed = 1
)
dev.off()

combined.seurat$nb <- factor(as.numeric(combined.seurat$Fate), levels = 1:15)
Idents(combined.seurat) <- "nb"
my.colors <- c("#1f77b4", "#5254a3", "#e377c2", "#ffbb78", "#e6550d",
               "#98df8a", "#f7b6d2", "#ff9896", "#9467bd", "#c5b0d5",
               "#8c564b", "#c49c94", "#bcbd22", "#dbdb8d", "#1bca8d")
names(my.colors) <- 15:1
pdf(file = file.path(plot.dir, "clusters.pdf"),
    width = 6 , height = 5)
DimPlot(combined.seurat,
        cols = my.colors,
        label = T,
        pt.size = 1.5, 
        label.size = 4,
        shuffle = T,
        seed = 1)
dev.off()
combined.seurat$nb_reverse_order <- factor(as.numeric(combined.seurat$Fate), levels = 15:1)
Idents(combined.seurat) <- "nb_reverse_order"

pdf(file = file.path(plot.dir, "dotplot.pdf"),
    width = 9.3 , height = 7.5)
DotPlot_scCustom(seurat_object = combined.seurat, features = c("T", "Sox2","Nkx1-2", "Fgf8",
                                                               "Foxa2", "Sox17","Msgn1",
                                                               "Tbx6", "Foxc1","Gata6","Tbx1",
                                                               "Irx3","Sox3",
                                                               "Cdx2","Cdx4","Fgf17","Epha5",
                                                               "Hes7", "Rspo3","Dll1" ,
                                                               "Mesp2","Ripply2","Cer1","Aldh1a2",
                                                               "Meox1","Tcf15","Pax3",
                                                               "Pax1", "Pax7","Igfbp5",
                                                               "Grsf1","Cdx1","Fgfbp3","Fgf3","Fgf4","Wnt6",
                                                               "Etv4",
                                                               "Klf2","Zfp42","Esrrb"),
                 colors_use = viridis_plasma_dark_high,
                 x_lab_rotate = T) +
  theme(legend.position = "bottom")
dev.off()

### number of cells per cluster
num_cells <- table(combined.seurat$Genotype, combined.seurat$Fate)
num_cells <- sweep(num_cells, 2, colSums(num_cells), "/")
num_cells <- as.data.frame.matrix(num_cells)
num_cells <- cbind(rownames(num_cells), data.frame(num_cells, row.names=NULL, check.names = FALSE))
colnames(num_cells)[1] <- "Genotype"
num_cells <- pivot_longer(num_cells,
                          cols = !matches("Genotype"), 
                          names_to = "Tissue",
                          values_to = "Percentage"
)
num_cells$Tissue <- factor(num_cells$Tissue, levels = levels(combined.seurat$Fate_reverse_order))
pdf(file = file.path(plot.dir, "pourcentage_cluster.pdf"),
    width = 8 , height = 5)
ggplot(data = num_cells, aes(x=Tissue, y=Percentage, fill = Genotype, factor())) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.6) +
  theme_classic() +
  scale_fill_manual(values=c("#e69a8d", 
                             "#5f4b8b")) +
  coord_flip()
dev.off()

num_cells_cluster_Time <- combined.seurat[[]] %>%
  group_by(nb, Time, Genotype) %>%
  summarize(n_cells = n()) %>%
  group_by(nb) %>%
  mutate(prop = n_cells / sum(n_cells))

num_cells_cluster_Time$Time <- factor(num_cells_cluster_Time$Time, levels = c("72h", "96h", "120h", "144h"))

ggplot(data = num_cells_cluster_Time, aes(x = prop, y = Genotype)) +
  geom_bar(stat = "identity", aes(fill = Genotype, alpha = Time)) +
  facet_grid(nb ~ ., switch = "y") +
  theme_classic() +
  scale_fill_manual(values=c("#e69a8d",
                             "#5f4b8b")) +
  scale_alpha_manual(values = c(0.5, .7, 1)) +
  scale_x_continuous("Fraction per cluster", expand = expansion(mult = c(0.001, .05))) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave(file = file.path(plot.dir, "pourcentage_cluster_timeinfo.pdf"),
    width = 3 , height = 5)

mapal <- viridis(64,direction = -1)

order.index <- order(combined.seurat$velocity_pseudotime, decreasing = F)
length(order.index)
# 17750
cell.order <- colnames(combined.seurat)[order.index]

# Get scale values for all genes:
combined.seurat <- ScaleData(combined.seurat, features = rownames(combined.seurat))
combined.seurat.NMP.Meso.pathway <- subset(combined.seurat,cells = WhichCells(combined.seurat, expression =
                                                                                (Time %in% c("120h") & (Fate %in% c("NMPs", "Anterior PSM", "Posterior PSM",
                                                                                             "Somitic Mesoderm", "Dermo/Sclero"))) |
                                                                                (Time %in% c("144h") & (Fate %in% c("Dermo/Sclero")))
                                                                              ))
table(combined.seurat.NMP.Meso.pathway$Fate)
table(combined.seurat.NMP.Meso.pathway$Genotype)
Idents(combined.seurat.NMP.Meso.pathway) <- "nb"
DimPlot_scCustom(
  combined.seurat.NMP.Meso.pathway,
  colors_use = my.colors[10:6],
  split.by = "Genotype",
  pt.size = 0.8, num_columns = 1) & xlim (c(-8, 5)) & NoAxes()
ggsave(file.path(plot.dir, "clusters_120h-144h_NMP_Meso.pdf"),
       width = 3.5 , height = 6)

FeaturePlot_scCustom(
  combined.seurat.NMP.Meso.pathway,
  features = "velocity_pseudotime",
  colors_use = cividis(64, direction = -1),
  split.by = "Genotype",
  alpha_exp = 1,
  na_cutoff = min(combined.seurat.NMP.Meso.pathway$velocity_pseudotime),
  pt.size = 0.8, num_columns = 1) &
  theme(legend.title = element_blank()) & xlim (c(-8, 5)) & NoAxes()
ggsave(file.path(plot.dir, "pseudotime_UMAP_120h-144h_NMP_Meso.pdf"),
       width = 3.5 , height = 6)

pdf(file = file.path(plot.dir, "doHeatmap_Genes_NMP_meso_only120h-144h.pdf"),
    width = 8 , height = 5)
DoHeatmap(combined.seurat.NMP.Meso.pathway,
          cells = intersect(cell.order,colnames(combined.seurat.NMP.Meso.pathway)),
          features = c("T", "Sox2","Nkx1-2","Cdx2","Cdx4", "Wnt3a", "Fgf8","Fgf17","Epha5",
                       "Crabp2","Hes7", "Rspo3","Tbx6","Msgn1","Dll1" , "Notch1",
                       "Mesp2","Ripply2","Aldh1a2","Lfng","Pcdh8", "Cer1","Snai1", "Twist1","Prrx2",
                       "Foxp1","Foxc1","Tcf15","Meox1","Foxd1","Pax3","Uncx",
                       "Igfbp5","Cxcl12","Meox2","Pax1","Pax7"),
          disp.min = -1,
          raster = F,
          group.by = "Genotype",
          slot = "scale.data") +
  scale_fill_gradientn(colours = rev(mapal))
dev.off()




# pseudotime ordred cell bar
mapal <- cividis(64,direction = -1)
combined.seurat.NMP.Meso.pathway[["pseudotime"]] <- CreateAssayObject(data = 
                                                                        t(x = FetchData(object = 
                                                                                          combined.seurat.NMP.Meso.pathway,
                                                                                        vars = 'velocity_pseudotime')))
pdf(file = file.path(plot.dir, "doHeatmap_pseudotime_meso.pdf"),
    width = 8 , height = 1)
DoHeatmap(combined.seurat.NMP.Meso.pathway,
          cells = intersect(cell.order,colnames(combined.seurat.NMP.Meso.pathway)),
          features = "velocity-pseudotime",
          disp.min = -1,
          raster = F,
          group.by = "Genotype",
          assay = 'pseudotime',
          slot = "data") +
  scale_fill_gradientn(colours = rev(mapal))
dev.off()

## gene expression for Cer1, Pcdh8, Pax1 and Pax7
for (gene in c("Cer1", "Pcdh8", "Pax1", "Pax7","Snai1")) {
  pdf(file = file.path(plot.dir, paste0("geneexpression_NMP.Meso.pathway_", gene, ".pdf")),
      width = 3.5 , height = 6)
  g <- FeaturePlot_scCustom(combined.seurat.NMP.Meso.pathway , features = gene,
                       split.by = "Genotype",
                       alpha_exp = 1,
                       pt.size = 0.8,
                       num_columns = 1) & xlim (c(-8, 5)) & NoAxes()
  print(g)
  dev.off()
}


## Hox expression at 72h, 96h, 120h and 144h
my.genes <- intersect(
  paste0("Hox", rep(letters[1:4], each = 13), rep(1:13, 4)),
  rownames(combined.seurat))
combined.seurat$genotype_time <- paste0(combined.seurat$Genotype,"_",combined.seurat$Time)

tmp <- PseudobulkExpression(
  object = combined.seurat,
  group.by = "genotype_time",
  # features = my.genes, # I want to normalize to all genes so I must not set features here
  method = "aggregate",
  layer = "count",
  return.seurat = TRUE
)
matrix.average.norm <- tmp[["RNA"]]$data[my.genes, ]

col_order <- c("wt-72h", "wt-96h", "wt-120h", "wt-144h",
               "BADC-72h", "BADC-96h", "BADC-120h","BADC-144h")
pheatmap(
  matrix.average.norm[,col_order],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  filename = file.path(plot.dir, "Heatmap_Hox.pdf")
)

# Specific genes
hoxa.b.c.d.genes <- read.delim(file.path(wd, "../../Code_figures/outputs/general/delBADC_genes.txt"), header = FALSE)$V1
chry.genes <- read.delim(file.path(wd, "../../Code_figures/outputs/general/chrYgenes.txt"), header = FALSE)$V1

Idents(combined.seurat) <- "Fate"

markers.ps <- FindMarkers(combined.seurat, ident.1 = "Primitive Streak", only.pos = TRUE)
selected.markers <- subset(markers.ps, p_val_adj < 0.01 & pct.1 > 0.2)
selected.markers <- selected.markers[order(-selected.markers$avg_log2FC), ]
selected.markers.genes <- sort(rownames(selected.markers)[1:50])
cat(selected.markers.genes)

# differential markers of primitive streak at 72h
combined.seurat.72h <- subset(combined.seurat,cells = WhichCells(combined.seurat, expression =
                                                                                   Time == "72h"))
combined.seurat.72h$genotype_fate <- paste0(combined.seurat.72h$Genotype,"_",combined.seurat.72h$Fate)
Idents(combined.seurat.72h) <- combined.seurat.72h$genotype_fate
markers.primitivestreak.wtvsBADC <- FindMarkers(combined.seurat.72h,
                                                ident.1 = "wt_Primitive Streak", ident.2 = "BADC_Primitive Streak",
                                                logfc.threshold = 0, min.pct = 0.01)
markers.primitivestreak.wtvsBADC$name <- rownames(markers.primitivestreak.wtvsBADC)
markers.primitivestreak.wtvsBADC$loci <- "other"
markers.primitivestreak.wtvsBADC$loci[rownames(markers.primitivestreak.wtvsBADC) %in% hoxa.b.c.d.genes] <- "HoxABCD"
markers.primitivestreak.wtvsBADC$loci[rownames(markers.primitivestreak.wtvsBADC) %in% chry.genes] <- "chrY"

write.table(markers.primitivestreak.wtvsBADC, file.path(plot.dir, "primitivestreak_72h_wtVSBADC.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
significant.genes <- subset(markers.primitivestreak.wtvsBADC, p_val_adj < 0.01 & (avg_log2FC > 1.5 | avg_log2FC < -1.5))

# table(significant.genes$name %in% Command(combined.seurat, "RunPCA.RNA")$features)
# FALSE  TRUE
# 202    42

pdf(file = file.path(plot.dir, "Volcano_primitivestraek.wtvsBADC_72h.pdf"),
    width = 10 , height = 10)
ggplot(data=markers.primitivestreak.wtvsBADC, aes(x=-avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(aes(color = loci)) + 
  theme_classic() +
  geom_text_repel(data=significant.genes,
                  aes(-avg_log2FC,-log10(p_val_adj),label=name)) +
  geom_point(data = subset(markers.primitivestreak.wtvsBADC, name %in% selected.markers.genes), shape = 1, col = "blue") +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red") +
  scale_color_manual("Gene location", values = c("other" = "black", "HoxABCD" = "pink", "chrY" = "cyan"))
dev.off()

# differential markers of each cluster at 120h
combined.seurat.120h <- subset(combined.seurat,cells = WhichCells(combined.seurat, expression =
                                                                   Time == "120h"))
combined.seurat.120h$genotype_fate <- paste0(combined.seurat.120h$Genotype,"_",combined.seurat.120h$Fate)
Idents(combined.seurat.120h) <- combined.seurat.120h$genotype_fate
for (fate in c("Posterior PSM", "Anterior PSM", "NMPs", "Somitic Mesoderm")) {
  markers.fate.wtvsBADC <- FindMarkers(
    combined.seurat.120h,
    ident.1 = paste0("wt_", fate), ident.2 = paste0("BADC_", fate),
    logfc.threshold = 0, min.pct = 0.01
  )
  markers.fate.wtvsBADC$name <- rownames(markers.fate.wtvsBADC)
  markers.fate.wtvsBADC$loci <- "other"
  markers.fate.wtvsBADC$loci[rownames(markers.fate.wtvsBADC) %in% hoxa.b.c.d.genes] <- "HoxABCD"
  markers.fate.wtvsBADC$loci[rownames(markers.fate.wtvsBADC) %in% chry.genes] <- "chrY"
  
  write.table(markers.fate.wtvsBADC, file.path(plot.dir, paste0(fate, "_120h_wtVSBADC.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  significant.genes = subset(markers.fate.wtvsBADC, p_val_adj < 0.01 & (avg_log2FC > 1.5 | avg_log2FC < -1.5))
  pdf(file = file.path(plot.dir, paste0("Volcano_", fate, ".wtvsBADC_120h.pdf")),
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
  
}
# expression of anterior primitive streak, endodrem and pharyngeal mesoderm
for (gene in c("T", "Gsc", "Foxa2", "Sox17","Tbx1","Gata6")) {
  pdf(file = file.path(plot.dir, paste0("geneexpression_APS_endomesoderm_", gene, ".pdf")),
      width = 3.5 , height = 6)
  g <- FeaturePlot_scCustom(combined.seurat , features = gene,
                            split.by = "Genotype",
                            alpha_exp = 1,
                            pt.size = 0.8,
                            num_columns = 1) & NoAxes()
  print(g)
  dev.off()
}

# expression of caudal epibast markers
for (gene in c("Cdx1", "Grsf1", "Etv4", "Fgf3","Fgf4")) {
  pdf(file = file.path(plot.dir, paste0("geneexpression_CaudalEpiblast_", gene, ".pdf")),
      width = 3.5 , height = 6)
  g <- FeaturePlot_scCustom(combined.seurat , features = gene,
                            split.by = "Genotype",
                            alpha_exp = 1,
                            pt.size = 0.5,
                            num_columns = 1) & NoAxes()
  print(g)
  dev.off()
}

# expression of caudal to epiblast markers in mutant gastruloid
combined.seurat.BADC <- subset(combined.seurat,cells = WhichCells(combined.seurat, expression =
                                                                    Genotype == "BADC"))

for (gene in c("T", "Nkx1-2", "Msgn1", "Tbx6")) {
  pdf(file = file.path(plot.dir, paste0("geneexpression_BADC_epiblastMesoderm_", gene, ".pdf")),
      width = 3.5 , height = 3)
  g <- FeaturePlot_scCustom(combined.seurat.BADC , features = gene,
                            split.by = "Genotype",
                            alpha_exp = 1,
                            pt.size = 0.5,
                            num_columns = 1) +
    xlim (c(-8, 5)) +
    ylim (c(-3, 11)) +
    NoAxes() +
    ggtitle(gene) +
    theme(legend.title = element_blank(), plot.title = element_text(face = "italic"))
  print(g)
  dev.off()
}

# expression of Cer1 in wt vs mutant gastruloid
pdf(file = file.path(plot.dir, paste0("geneexpression_Cer1.pdf")),
    width = 5 , height = 8)
g <- FeaturePlot_scCustom(combined.seurat , features = "Cer1",
                          split.by = "Genotype",
                          alpha_exp = 1,
                          pt.size = 0.5,
                          num_columns = 1)
print(g)
dev.off()

# expression of Wnt, Fgf and Notch markers
for (gene in c("Wnt3a")) {
  pdf(file = file.path(plot.dir, paste0("geneexpression_NMP.Meso.pathway_", gene, ".pdf")),
      width = 7 , height = 12)
  g <- FeaturePlot_scCustom(combined.seurat , features = gene,
                            split.by = "Genotype",
                            alpha_exp = 1,
                            pt.size = 1.2,
                            num_columns = 1) & NoAxes()
  print(g)
  dev.off()
}

