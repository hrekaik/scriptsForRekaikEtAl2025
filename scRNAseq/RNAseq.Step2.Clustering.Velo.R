# Install required packages
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("SeuratWrappers")
safelyLoadAPackageInCRANorBioconductor("sceasy")
safelyLoadAPackageInCRANorBioconductor("reticulate")
use_python("/usr/bin/python3")


# introduce variable that will be used
wd <- "/Papers/DelHox/toGEO/scRNAseq"
plot.dir <- "../../Code_figures/outputs/scRNAseq/"
setwd(wd)
dir.create(plot.dir, showWarnings = FALSE)
# all files are relative to wd
# rds files will be put in a RDS folder in the working directory
# rds are only created once, if pipeline changes, delete manually all rds
RDSfolder <- "RDS"

#here put name of RDS as given in step 2

nameRDS <- "combined_AlltimewtvsBADC_96f83823.RDS"
combined.seurat <- readRDS(file.path(RDSfolder, paste0(nameRDS)))



ElbowPlot(combined.seurat, ndims = 50)
seed <- 5
combined.seurat <- RunUMAP(combined.seurat, reduction = "pca",
                           n.components = 2L,  
                           dims = 1:30, seed.use = seed)
combined.seurat <- FindNeighbors(combined.seurat, reduction = "pca",
                                 dims = 1:30)
# Default resolution
combined.seurat <- FindClusters(combined.seurat, resolution = 0.80)

combined.seurat <- FindClusters(combined.seurat, resolution = 1.5)

combined.seurat <- JoinLayers(combined.seurat)
DimPlot(combined.seurat, group.by = "Genotype") +
  DimPlot(combined.seurat, group.by = "Time") +
  DimPlot(combined.seurat, group.by = "RNA_snn_res.0.8", label = TRUE) +
  DimPlot(combined.seurat, group.by = "RNA_snn_res.1.5", label = TRUE)
ggsave(file.path(plot.dir, "clusters_1.5.pdf"), width = 10, height = 10)

all.cl <- FindAllMarkers(combined.seurat, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.2)
write.csv(all.cl, file.path(plot.dir, "clusters_1.5_markers_stringent.csv"))
top2.markers <- all.cl %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)
DotPlot(combined.seurat, features = unique(top2.markers$gene)) + RotatedAxis()
ggsave(file.path(plot.dir, "clusters_1.5_dotplot2.pdf"), width = 15, height = 10)

#### clusters name
### Change clusters name and order
Idents(combined.seurat) <- combined.seurat$seurat_clusters


new.cluster.name <-  c(
  "0" = "Nascent Mesoderm",
  "1" = "Primitive Streak",
  "2" = "Somitic Mesoderm",
  "3" = "Primitive Streak",
  "4" = "Dermo/Sclero",
  "5" = "NMPs",
  "6" = "Caudal Epiblast",
  "7" = "Neural Tube",
  "8" = "Caudal Mesoderm",
  "9" = "Dermo/Sclero",
  "10" = "Anterior PSM",
  "11" = "Mixed Mesoderm",
  "12" = "Neural Tube",
  "13" = "Posterior PSM",
  "14" = "Definitive Endoderm",
  "15" = "PGC",
  "16" = "Mixed Mesoderm",
  "17" = "Somitic Mesoderm",
  "18" = "Unknown",
  "19" = "Endothelial"
)


combined.seurat <- RenameIdents(combined.seurat, new.cluster.name)
combined.seurat[["Fate"]] <- Idents(combined.seurat)

DimPlot(combined.seurat, label = TRUE)

levels(x=combined.seurat) <- c("Unknown","Endothelial","PGC",
                               "Caudal Mesoderm","Caudal Epiblast",
                               "Dermo/Sclero","Somitic Mesoderm","Anterior PSM","Posterior PSM","NMPs", "Neural Tube",
                               "Mixed Mesoderm","Nascent Mesoderm","Definitive Endoderm","Primitive Streak")
combined.seurat$Fate_reverse_order <- factor(as.character(combined.seurat$Fate), levels = levels(combined.seurat))
combined.seurat$Fate <- factor(as.character(combined.seurat$Fate), levels = rev(levels(combined.seurat)))
DimPlot(combined.seurat, group.by = "Fate", label = TRUE)
saveRDS(combined.seurat, file.path(RDSfolder, nameRDS))
### Velocyto

combined.seurat.novelo <- combined.seurat
# The barcodes in the loom file are like: "72h_WT_rep1:AACCACAAGCAATTCCx" "72h_WT_rep1:AACAACCAGTACTCGTx" "72h_WT_rep1:AACAAGAAGTGGTTAAx"
# While in Seurat they are like: "72h.WT.rep1_AAACCCAAGCTGGCCT" "72h.WT.rep1_AAACCCACATGACTCA" "72h.WT.rep1_AAACGAAAGCCACAAG"
# Generate a named vector for conversion:
combined.seurat.barcodes <- sapply(strsplit(colnames(combined.seurat.novelo), "_"), tail, 1)
combined.seurat.barcodes.conversion <- paste0(combined.seurat.novelo$RDS, "_", combined.seurat.barcodes)
names(combined.seurat.barcodes.conversion) <- paste0(combined.seurat.novelo$orig.ident, ":",
                                                     combined.seurat.barcodes, "x")
LoomFolder <- "Raw/loom"
# Make vector with loom paths
loom.files <- file.path(LoomFolder, paste0(unique(combined.seurat.novelo$orig.ident), ".loom"))
names(loom.files) <- unique(combined.seurat.novelo$orig.ident)
# Read them and replace name with
my.function.readandreplacename.velocity <- function(names.loom, loom) {
  loom.matrix <- ReadVelocity(loom)
  cell.names <- colnames(loom.matrix$spliced)
  # Initial cell names are:
  # old_name:barcodex
  new.cell.names <- paste0(names.loom, ":", (sapply(strsplit(cell.names, ":"), tail, 1)))
  # Replace old_name by the name of the cellplex
  colnames(loom.matrix$spliced) <- new.cell.names
  stopifnot(all(colnames(loom.matrix$unspliced) == cell.names))
  colnames(loom.matrix$unspliced) <- new.cell.names
  stopifnot(all(colnames(loom.matrix$ambiguous) == cell.names))
  colnames(loom.matrix$ambiguous) <- new.cell.names
  names(loom.matrix) 
  return(loom.matrix)
}
matrices <- lapply(seq_along(loom.files),
                   function(i) my.function.readandreplacename.velocity(names(loom.files)[i],unname(loom.files)[i]))

# remove duplicated gene names
my.function.remove.duplicates <- function(mtx) {
  original.names <- rownames(mtx$spliced)
  new.names <- make.unique(original.names)
  rownames(mtx$spliced) <- new.names
  stopifnot(all(rownames(mtx$unspliced) == original.names))
  rownames(mtx$unspliced) <- new.names
  stopifnot(all(rownames(mtx$ambiguous) == original.names))
  rownames(mtx$ambiguous) <- new.names
  return(mtx)
}
matrices_clean <- lapply(matrices, my.function.remove.duplicates)

rm(matrices)

# Convert to Seurat to more easily rename cells and genes and transfer assays:
seurats <- lapply(matrices_clean, as.Seurat)
names(seurats) <- names(loom.files)

rm(matrices_clean)

# Combine Seurat objects keeping the same order:
Combined.seurat.velocity <- merge(seurats[[as.character(unique(combined.seurat.novelo$orig.ident))[1]]],
                                  y = subsetByNamesOrIndices(seurats,
                                                             as.character(unique(combined.seurat.novelo$orig.ident))[2:length(as.character(unique(combined.seurat.novelo$orig.ident)))]),
                                  project = "aggregate.velocity")
rm(seurats)
Combined.seurat.velocity <- JoinLayers(Combined.seurat.velocity)
# Restrict to cells which are on the combined.seurat object:
Combined.seurat.velocity <- subset(Combined.seurat.velocity, cells = names(combined.seurat.barcodes.conversion))
# Rename cells
Combined.seurat.velocity <- RenameCells(Combined.seurat.velocity, new.names = unname(combined.seurat.barcodes.conversion[colnames(Combined.seurat.velocity)]))
# Convert v5 to v3
Combined.seurat.velocity[["spliced"]] <- CreateAssayObject(counts = Combined.seurat.velocity@assays$spliced$counts)

# Transfer assays
combined.seurat.novelo[["spliced"]] <- Combined.seurat.velocity[["spliced"]]
combined.seurat.novelo[["unspliced"]] <- Combined.seurat.velocity[["unspliced"]]


temp <- DietSeurat(combined.seurat.novelo,
                   counts = TRUE,
                   data = TRUE,
                   scale.data = FALSE,
                   features = NULL,
                   assays = NULL,
                   dimreducs = Reductions(combined.seurat),
                   misc = TRUE)


# Restrict to genes in RNA:
rna.genes <- rownames(temp[["RNA"]])
temp[["spliced"]] <- subset(temp[["spliced"]], features = rna.genes)
temp[["unspliced"]] <- subset(temp[["unspliced"]], features = rna.genes)
# Convert v5 to v3
temp[["RNA"]] <- as(temp[["RNA"]], Class = "Assay")


# export to annadata
sceasy::convertFormat(temp, from="seurat", to="anndata",
                      outFile=file.path("Hocine_prevelo_diet.h5ad"))
write.csv(x = t(as.matrix(temp@assays$spliced@counts)), file = file.path('Hocine_prevelo_diet_spliced.csv'))
write.csv(x = t(as.matrix(temp@assays$unspliced@counts)), file = file.path('Hocine_prevelo_diet_unspliced.csv'))


system(paste("/usr/bin/python3 ../../Code_figures/scRNAseq/scRNA-seq_scvelo.py", wd))
system(paste("/usr/bin/python3 ../../Code_figures/scRNAseq/scRNA-seq_plot_scvelo.py", wd, plot.dir))


### open h5ad file and use velocity data
seurat.velocity <- sceasy::convertFormat("Hocine_velo_diet.h5ad",
                                         from="anndata", to="seurat")
all(colnames(seurat.velocity) == colnames(combined.seurat))
combined.seurat$velocity_pseudotime <- seurat.velocity$velocity_pseudotime
combined.seurat$latent_time <- seurat.velocity$latent_time

# save RDS file
saveRDS(combined.seurat, gsub(".RDS", ".analyzed.scVelo.RDS", nameRDS))
