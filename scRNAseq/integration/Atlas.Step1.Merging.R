options(timeout = 600)

library(usefulLDfunctions)

safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("ggplot2")

wd <- "/Papers/DelHox/toGEO/scRNAseq"
dir.with.mouse.RDS <- "../../RDS_mouse_TOME"
dir.github <- "../../Code_figures"
hoxa.b.c.d.genes <- read.delim(file.path(wd, dir.github, "outputs/general/delBADC_genes.txt"), header = FALSE)$V1

setwd(wd)

nameRDS <- "combined_AlltimewtvsBADC_96f83823.analyzed.scVelo.RDS"
rds.atlas <- file.path(dir.with.mouse.RDS, "seurat.atlas.combined.E6.5-9.5.rds")
if (!file.exists(rds.atlas)) {
  ###### Downlad and read available Atals RDS files (https://tome.gs.washington.edu)
  urlfiles <- c(
    "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/tome_summary_data/mm/seurat_object_E6.5.rds", # Mohammed et al., Cheng et al., Pijuan-Sala et al.
    "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/tome_summary_data/mm/seurat_object_E6.75.rds", # Pijuan-Sala et al.
    "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/tome_summary_data/mm/seurat_object_E7.0.rds", # Pijuan-Sala et al.
    "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/tome_summary_data/mm/seurat_object_E7.25.rds", # Pijuan-Sala et al.
    "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/tome_summary_data/mm/seurat_object_E7.5.rds", # Pijuan-Sala et al.
    "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/tome_summary_data/mm/seurat_object_E7.75.rds", # Pijuan-Sala et al.
    "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/tome_summary_data/mm/seurat_object_E8.rds", # Pijuan-Sala et al.
    "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/tome_summary_data/mm/seurat_object_E8.25.rds", # Pijuan-Sala et al.
    "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/tome_summary_data/mm/seurat_object_E8.5a.rds", # Pijuan-Sala et al.
    "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/tome_summary_data/mm/seurat_object_E9.5.rds" # Cao et al.
  )
  dir.create(dir.with.mouse.RDS, recursive = TRUE, showWarnings = FALSE)
  lapply(urlfiles, function(url) {
    output.file <- file.path(dir.with.mouse.RDS, basename(url))
    if (!file.exists(output.file)) {
      download.file(url, output.file)
    }
  })

  my.stages <- c("E6.5", "E6.75", "E7.0", "E7.25", "E7.5", "E7.75", "E8", "E8.25", "E8.5a", "E9.5")

  ##### combine all atlas stage in one seurat object and keep only maximum 20k cells
  seurat.atlas.subset.combined <- list()
  for (my.stage in my.stages) {
    cat(paste0("loading ", my.stage, "\n"))
    seurat.atlas <- readRDS(file.path(dir.with.mouse.RDS, paste0("seurat_object_", my.stage, ".rds")))
    Idents(seurat.atlas) <- "day"
    cat(paste0("subset ", my.stage, "\n"))
    seurat.atlas.subset <- subset(seurat.atlas, downsample = 20000)
    seurat.atlas.subset@meta.data[["old.ident"]] <- "Mouse"
    seurat.atlas.subset.combined[[my.stage]] <- seurat.atlas.subset
  }
  # Free memory
  rm(list = c("seurat.atlas", "seurat.atlas.subset"))
  gc()

  lapply(seurat.atlas.subset.combined, dim)
  # $E6.5
  # [1] 31435  4444
  #
  # $E6.75
  # [1] 29452  2075
  #
  # $E7.0
  # [1] 29452 14749
  #
  # $E7.25
  # [1] 29452 13537
  #
  # $E7.5
  # [1] 29452 10994
  #
  # $E7.75
  # [1] 29452 14493
  #
  # $E8
  # [1] 29452 16681
  #
  # $E8.25
  # [1] 29452 15935
  #
  # $E8.5a
  # [1] 29452 16909
  #
  # $E9.5
  # [1] 24552 20000

  # Merge
  seurat.atlas.combined <- merge(
    x = seurat.atlas.subset.combined[[1]],
    y = seurat.atlas.subset.combined[2:length(seurat.atlas.subset.combined)],
    project = "aggregate.all"
  )

  # Free memory space
  rm(seurat.atlas.subset.combined)
  gc()
  ### Replace gene ID with gene symbol and tag unidentified genes
  conversion.table <- read.delim(file.path(dir.github, "outputs", "scRNAseq", "integration", "gene_id_name_conversion.txt"))
  all(rownames(seurat.atlas.combined) %in% conversion.table$gene_ID)
  conversion.name <- conversion.table$name
  names(conversion.name) <- conversion.table$gene_ID

  seurat.atlas.combined@assays$RNA@counts@Dimnames[[1]] <- conversion.name[seurat.atlas.combined@assays$RNA@counts@Dimnames[[1]]]
  seurat.atlas.combined@assays$RNA@data@Dimnames[[1]] <- conversion.name[seurat.atlas.combined@assays$RNA@data@Dimnames[[1]]]
  rownames(seurat.atlas.combined[["RNA"]]@meta.features) <- conversion.name[rownames(seurat.atlas.combined[["RNA"]]@meta.features)]
  seurat.atlas.combined <- NormalizeData(seurat.atlas.combined, normalization.method = "LogNormalize")
  saveRDS(seurat.atlas.combined, file = rds.atlas)
} else {
  seurat.atlas.combined <- readRDS(rds.atlas)
}

# Get the gastruloid seurat
seurat.gastruloid <- readRDS(nameRDS)
seurat.gastruloid@meta.data[["old.ident"]] <- "Gastruloid"
seurat.gastruloid$origin <- "Gastruloid"
  
# Run integration using https://satijalab.org/seurat/articles/seurat5_integration
# First get a single seurat object
seurat.both <- merge(seurat.atlas.combined, seurat.gastruloid)
rm(list = c("seurat.atlas.combined", "seurat.gastruloid"))
gc()
# Split layers by study
seurat.both[["RNA"]] <- JoinLayers(seurat.both[["RNA"]])
gc()
seurat.both$study <- "our"
seurat.both$study[seurat.both$old.ident == "Mouse"] <- "TOME_Pijuan-Sala"
seurat.both$study[seurat.both$group %in% "352"] <- "TOME_Cheng"
seurat.both$study[seurat.both$group %in% "351"] <- "TOME_Mohammed"
seurat.both$study[seurat.both$group %in% "E9.5"] <- "TOME_E9.5"
table(seurat.both$study, seurat.both$old.ident)
table(seurat.both$study, seurat.both$day, exclude = NULL)
seurat.both[["RNA"]] <-  split(seurat.both[["RNA"]], f = seurat.both$study)
gc()
seurat.both <- NormalizeData(seurat.both)
seurat.both <- FindVariableFeatures(seurat.both, selection.method = "vst", nfeatures = 2500)
VariableFeatures(seurat.both) <- setdiff(VariableFeatures(seurat.both), hoxa.b.c.d.genes)
gc()
# Should we use Cell cycle?? It seems that the cells do not cluster by Phase
seurat.both <- ScaleData(seurat.both)
seurat.both <- RunPCA(seurat.both)
gc()
# Use CCA:
# Impose the sample.tree
# Samples are TOME_Mohammed, TOME_Cheng, TOME_Pijuan-Sala, TOME_E9.5, our
sample.tree <- matrix(
  c(
    -3, -1, # Pijuan-Sala with Mohammed
    1, -2, # Add Cheng
    2, -4, # Add E9.5
    3, -5 # Add gastruloids
  ),
  byrow = TRUE, ncol = 2
)
seurat.both <- IntegrateLayers(
  object = seurat.both, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca.study",
  features = VariableFeatures(seurat.both),
  sample.tree = sample.tree
)
# Merging dataset 1 into 3
# Merging dataset 2 into 3 1
# Merging dataset 4 into 3 1 2
# Merging dataset 5 into 3 1 2 4
seurat.both <- RunUMAP(seurat.both, reduction = "integrated.cca.study", dims = 1:30, reduction.name = "umap.cca.per.study")
DimPlot(seurat.both, reduction = "umap.cca.per.study", group.by = c("day", "study", "cell_type"))
saveRDS(seurat.both, file = file.path(dir.with.mouse.RDS, "seurat.atlas.combined.E6.5-9.5.gastruloids.integrated.study.rds"))
