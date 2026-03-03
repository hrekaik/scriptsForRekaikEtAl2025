if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")


library(devtools)
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("Signac")
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
safelyLoadAPackageInCRANorBioconductor("biovizBase") # To build annotations
# This package has the p6 version which we did not used
# safelyLoadAPackageInCRANorBioconductor("BSgenome.Mmusculus.UCSC.mm10")
safelyLoadAPackageInCRANorBioconductor("Rsamtools")
safelyLoadAPackageInCRANorBioconductor("BSgenome")
safelyLoadAPackageInCRANorBioconductor("ggplot2") # To save diagnostic plots
safelyLoadAPackageInCRANorBioconductor("stringr")
safelyLoadAPackageInCRANorBioconductor("glmGamPoi") # To run SCTransform

wd <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures"
setwd(wd)
plot.dir <- "outputs/Multiome/120h"
hoxa.b.c.d.genes.file <- "outputs/general/delBADC_genes.txt"


counts_wt <- Read10X("../toGEO/Multiome/120h/wt_120h_GEX")
fragpath_wt <- "../toGEO/Multiome/120h/wt_120h_ATAC.bgzip"

counts_BADC <- Read10X("../toGEO/Multiome/120h/BADC_120h_GEX")
fragpath_BADC <- "../toGEO/Multiome/120h/BADC_120h_ATAC.bgzip"

# get unique peaks from MACS2
unique_peaks <- import("../toGEO/Multiome/120h/120h_merged_peaks.bed")
genome.fasta <- "/Users/hocine.rekaik/Desktop/EXP/Gastruloids/scMultiome/21012025/data/mm10_UCSC.fa"
if (! file.exists(paste0(genome.fasta, ".fai"))) {
  indexFa(genome.fasta)
}
genome <- FaFile(genome.fasta)
genome.seqinfo <- seqinfo(genome)

# Make sure all peaks end at the end of chromosome size
end(unique_peaks) <-
  apply(
    cbind(
      end(unique_peaks),
      seqlengths(genome.seqinfo)[as.vector(seqnames(unique_peaks))]
    ),
    1,
    min
  )

# get gene annotations for mm10
# This file comes from https://zenodo.org/records/10079673
annotation <- import.gff("data/mm10_allGastruloids_min10_extended.gtf")
# See https://github.com/stuart-lab/signac/issues/1793
annotation$tx_id <- annotation$transcript_id

# GEX
## wt
seurat_wt <- CreateSeuratObject(
  counts = counts_wt,
  assay = "RNA"
)
## BADC
seurat_BADC <- CreateSeuratObject(
  counts = counts_BADC,
  assay = "RNA"
)

# ATAC
## wt
fragment_wt <- CreateFragmentObject(fragpath_wt)
peaks_ATAC_wt <- FeatureMatrix(
  fragments = fragment_wt,
  features = unique_peaks,
  cells = colnames(seurat_wt)
)

## BADC
fragment_BADC <- CreateFragmentObject(fragpath_BADC)
peaks_ATAC_BADC <- FeatureMatrix(
  fragments = fragment_BADC,
  features = unique_peaks,
  cells = colnames(seurat_BADC)
)
# Warning message:
# In SingleFeatureMatrix(fragment = fragments[[x]], features = features,  :
#   2 features are on seqnames not present in the fragment file. These will be removed.

# Check which they are:
missing.peaks <- setdiff(rownames(peaks_ATAC_wt), rownames(peaks_ATAC_BADC))
print(missing.peaks)
# 2 from chr4_JH584295_random
# Add them back
peaks_ATAC_BADC <- Signac:::AddMissing(
  peaks_ATAC_BADC,
  cells = colnames(seurat_BADC),
  features = Signac:::GRangesToString(unique_peaks)
)

seurat_wt[["ATAC"]] <- CreateChromatinAssay(
  counts = peaks_ATAC_wt,
  sep = c(":", "-"),
  fragments = fragment_wt,
  annotation = annotation,
  genome = genome.seqinfo
)
seurat_wt <- AddMetaData(seurat_wt, metadata = "wt", col.name = "Genotype")

seurat_BADC[["ATAC"]] <- CreateChromatinAssay(
  counts = peaks_ATAC_BADC,
  sep = c(":", "-"),
  fragments = fragment_BADC,
  annotation = annotation,
  genome = genome.seqinfo
)
seurat_BADC <- AddMetaData(seurat_BADC, metadata = "BADC", col.name = "Genotype")

#### combine
combined.seurat <- merge(seurat_wt, y = seurat_BADC, add.cell.ids = c("wt", "BADC"), project = "Multiome")

DefaultAssay(combined.seurat) <- "ATAC"

combined.seurat <- NucleosomeSignal(combined.seurat)

combined.seurat <- TSSEnrichment(combined.seurat)

DensityScatter(combined.seurat, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE, )

Idents(combined.seurat) <- combined.seurat$Genotype
VlnPlot(
  object = combined.seurat,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
ggsave(file.path(plot.dir, "ATAC_before_filtering.png"), width = 8, height = 4)
# filter out low quality cells
print(combined.seurat)
combined.seurat <- subset(
  x = combined.seurat,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1800 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
print(combined.seurat)
VlnPlot(
  object = combined.seurat,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
ggsave(file.path(plot.dir, "ATAC_after_filtering.png"), width = 8, height = 4)


DefaultAssay(combined.seurat) <- "RNA"
combined.seurat <- JoinLayers(combined.seurat)

s.genes <- str_to_title(cc.genes.updated.2019$s.genes)
g2m.genes <- str_to_title(cc.genes.updated.2019$g2m.genes)
combined.seurat <- NormalizeData(combined.seurat)
combined.seurat <- CellCycleScoring(combined.seurat,
                                    s.features = s.genes,
                                    g2m.features = g2m.genes,
                                    verbose = FALSE)

combined.seurat <- PercentageFeatureSet(combined.seurat, pattern = "^mt-", col.name = "percent.mt")
options(future.globals.maxSize = 600 * 1024^2)
combined.seurat <- SCTransform(combined.seurat,vars.to.regress = c("percent.mt", "S.Score",
                                                                   "G2M.Score"))
hoxa.b.c.d.genes <- read.delim(hoxa.b.c.d.genes.file, header = FALSE)$V1
manual.removal.genes <- hoxa.b.c.d.genes

# Reduce batch effect
average.expression <- rowMeans(AverageExpression(combined.seurat, assay = "SCT", layer = "counts")[[1]])
genes.to.exclude <- names(average.expression[average.expression < quantile(average.expression,
                                                                           0.05) | average.expression > quantile(average.expression,
                                                                                                                 0.8)])
used.var.gene <- setdiff(VariableFeatures(combined.seurat, assay = "SCT"),
                         c(genes.to.exclude, manual.removal.genes))



combined.seurat <- RunPCA(combined.seurat, features = used.var.gene)

DefaultAssay(combined.seurat) <- "ATAC"
combined.seurat <- FindTopFeatures(combined.seurat, min.cutoff = 5)
combined.seurat <- RunTFIDF(combined.seurat)
# remove Hox clusters feeatures from SVD computing
ATAC.peaks <- granges(combined.seurat)
Hox.coordinates <- GRanges(
  seqnames = c('chr6', "chr11", "chr15", "chr2"),
  ranges = IRanges(start = c(52153796, 96191745, 102916767, 74667800),
                   end = c(52263152, 96372449, 103042756, 74767365)),
  strand = c('*', '*', '*', '*')
)
Hox.features.rows <- as.matrix(GenomicRanges::findOverlaps(ATAC.peaks,Hox.coordinates))[,1]
Hox.features <- rownames(combined.seurat)[Hox.features.rows]
used.var.peaks <- setdiff(
  VariableFeatures(combined.seurat, assay = "ATAC"),
  Hox.features
)
print(length(used.var.peaks))
print(dim(combined.seurat))
combined.seurat <- RunSVD(combined.seurat, features = used.var.peaks)

# build a joint neighbor graph using both assays
combined.seurat <- FindMultiModalNeighbors(
  object = combined.seurat,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:40),
  modality.weight.name = c("SCT.weight", "ATAC.weight"),
  verbose = TRUE
)




saveRDS(combined.seurat, "../toGEO/Multiome/120h/combined_Multiome_120h.RDS")

