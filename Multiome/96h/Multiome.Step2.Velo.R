library(devtools)
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("SeuratWrappers")
safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("sceasy")
safelyLoadAPackageInCRANorBioconductor("reticulate")
use_python("/usr/bin/python3")

wd <- "/scratch/ldelisle/rstudio_test/"
setwd(wd)

multiome.seurat <- readRDS("combined_Multiome.RDS")
multiome.seurat$Genotype <- factor(x = multiome.seurat$Genotype, levels = c("wt","BADC"))

# The barcodes in the loom file are like: "wt:AAAGCACCACACCAACx" "BADC:AAGACATAGTCACTAAx"
# While in Seurat they are like: "wt_AAAGCACCACACCAAC" "BADC_AAGACATAGTCACTAA"
# Generate a named vector for conversion:
combined.seurat.barcodes <- sapply(strsplit(colnames(multiome.seurat), "_"), tail, 1)
combined.seurat.barcodes.conversion <- colnames(multiome.seurat)
names(combined.seurat.barcodes.conversion) <- paste0(as.character(multiome.seurat$Genotype), ":",
                                                     combined.seurat.barcodes, "x")
LoomFolder <- "inputs_galaxy/Loom"
# Make vector with loom paths
loom.files <- file.path(LoomFolder, paste0(unique(multiome.seurat$Genotype), ".loom"))
names(loom.files) <- unique(multiome.seurat$Genotype)
matrices <- lapply(loom.files, ReadVelocity)
# Restrict to genes in RNA and reorder them:
rna.genes <- rownames(multiome.seurat[["RNA"]])
matrices.reordered <- lapply(matrices, function(m) {
  lapply(m, function(mat) {
    rownames(mat) <- make.unique(rownames(mat))
    mat[rna.genes, ]
  })
})
# seurats <- lapply(matrices.reordered, as.Seurat)
# The command above is failing because of the duplicated rownames
options(Seurat.object.assay.version = "v3")
seurats <- lapply(matrices.reordered, as.Seurat)
# Combine Seurat objects keeping the same order:
Combined.seurat.velocity <- merge(seurats[[as.character(unique(multiome.seurat$Genotype))[1]]],
                                  y = subsetByNamesOrIndices(seurats,
                                                             as.character(unique(multiome.seurat$Genotype))[2:length(as.character(unique(multiome.seurat$Genotype)))]),
                                  project = "aggregate.velocity")

# Restrict to cells which are on the subset.seurat object:
Combined.seurat.velocity <- subset(Combined.seurat.velocity, cells = names(combined.seurat.barcodes.conversion))

# Rename cells
Combined.seurat.velocity <- RenameCells(Combined.seurat.velocity, new.names = unname(combined.seurat.barcodes.conversion[colnames(Combined.seurat.velocity)]))
# Transfer assays
multiome.seurat[["spliced"]] <- Combined.seurat.velocity[["spliced"]]
multiome.seurat[["unspliced"]] <- Combined.seurat.velocity[["unspliced"]]

DefaultAssay(multiome.seurat) <- "RNA"
diet.multiome.seurat <- DietSeurat(multiome.seurat,
                                   counts = TRUE,
                                   data = TRUE,
                                   scale.data = FALSE,
                                   features = NULL,
                                   assays = c("RNA", "spliced", "unspliced"),
                                   dimreducs = Reductions(multiome.seurat),
                                   misc = TRUE)
# Convert v5 to v3:
diet.multiome.seurat[["RNA"]] <- CreateAssayObject(counts = GetAssayData(diet.multiome.seurat[["RNA"]], layer = "counts"))
SeuratDisk::SaveH5Seurat(
  diet.multiome.seurat,
  filename = "Multiome_prevelo.h5Seurat",
  overwrite = TRUE
)
SeuratDisk::Convert("Multiome_prevelo.h5Seurat", dest = "h5ad", overwrite = TRUE)

# Also export a version with less cells:
diet.multiome.seurat.caudal <- subset(diet.multiome.seurat, seurat_clusters %in% 0:3)
table(diet.multiome.seurat.caudal$seurat_clusters)
#    0    1    2    3    4    5    6 
# 2762 2393 1104  993    0    0    0 
SeuratDisk::SaveH5Seurat(
  diet.multiome.seurat.caudal,
  filename = "Multiome_caudal_prevelo.h5Seurat",
  overwrite = TRUE
)
SeuratDisk::Convert("Multiome_caudal_prevelo.h5Seurat", dest = "h5ad", overwrite = TRUE)

system(paste("/usr/bin/python3 scriptsForRekaikEtAl2025/Multiome/Multiome_scvelo.py", wd))

### open h5ad file and use velocity data
multiome.velo <- sceasy::convertFormat("Multiome_velo.h5ad",
                                       from = "anndata", to = "seurat")
multiome.seurat$velocity_pseudotime <- NA
multiome.seurat$velocity_pseudotime[match(colnames(multiome.velo), colnames(multiome.seurat))] <- multiome.velo$velocity_pseudotime
multiome.seurat$latent_time <- NA
multiome.seurat$latent_time[match(colnames(multiome.velo), colnames(multiome.seurat))] <- multiome.velo$latent_time

multiome.caudal.velo <- sceasy::convertFormat("Multiome_caudal_velo.h5ad",
                                              from = "anndata", to = "seurat")
multiome.seurat$caudal_velocity_pseudotime <- NA
multiome.seurat$caudal_velocity_pseudotime[match(colnames(multiome.caudal.velo), colnames(multiome.seurat))] <- multiome.caudal.velo$velocity_pseudotime
multiome.seurat$caudal_latent_time <- NA
multiome.seurat$caudal_latent_time[match(colnames(multiome.caudal.velo), colnames(multiome.seurat))] <- multiome.caudal.velo$latent_time

# save RDS file
saveRDS(multiome.seurat, "combined_Multiome_Velo.RDS")
