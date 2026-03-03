library(pheatmap)

gitHubDir <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures/RNAseq/RNAseq_ActA/ActA_wtvsMutant"
output_dir <- file.path(gitHubDir,"../../../outputs/RNAseq/RNAseq_ActA/ActA_wtvsMutant")

fpkms <- read.delim(file.path(gitHubDir,"../../../outputs/RNAseq/RNAseq_ActA/ActA_wtvsMutant/AllCufflinks_Simplified.txt"))

samples.plans <- read.delim(samplesPlan<-file.path(gitHubDir,"samplesPlan_ActA_wtvsBADC.txt"))





gene_sets <- list(
  Primitive_streak = c("T","Fgf8","Wnt3","Eomes"),
  Anterior_primitive_streak = c("Lhx1","Gsc","Cer1","Otx2"),
  Cardio_Pharyngeal_Mesoderm = c("Isl1","Mef2c","Nkx2-5","Myl7","Tcf21"),
  Def_Endoderm = c("Foxa2", "Gata6","Gata4","Sox17"),
  Cranial_mesoderm = c("Tbx1","Alx1","Ebf3"),
  Notochord = c("Noto","Chrd","Shh","Nog"),
  Paraxial_mesoderm = c("Msgn1","Tbx6" , "Dll1","Foxc2"),
  HoxLess_specific = c("Hoxb1", "Cdx1", "Grsf1")
)



rownames(samples.plans) <- samples.plans$sample
samples.plans$Replicate <- factor(samples.plans$Replicate)
samples.plans$Genotype <- samples.plans$Line



for (set_name in names(gene_sets)) {
  
  genes <- gene_sets[[set_name]]
  
  # Subset expression matrix
  mat <- fpkms[match(genes, fpkms$gene_short_name),
               paste0("FPKM_", samples.plans$sample)]
  
  # Remove missing genes
  keep <- !is.na(rownames(mat))
  mat <- mat[keep, , drop = FALSE]
  
  colnames(mat) <- samples.plans$sample
  rownames(mat) <- genes[keep]
  
  # Skip if no genes found
  if (nrow(mat) == 0) {
    warning(paste("No genes found for", set_name))
    next
  }
  
  # Output file
  out_file <- file.path(
    output_dir,
    paste0("heatmap_", set_name, ".pdf")
  )
  
  # Plot
  pheatmap(
    log2(1 + mat),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    main = set_name,
    annotation_col = samples.plans[, c("Replicate", "Genotype", "Tissue")],
    annotation_colors = list(
      "Genotype" = c("BADC" = "#e69a8d", "wt" = "#5f4b8b"),
      "Replicate" = c("1" = "white", "2" = "grey70"),
      "Tissue" = c("96h" = "orange", "120h" = "#606c38")
    ),
    file = out_file,
    width = 6,
    height = 1+ (0.3 * nrow(mat)) # auto-scale height
  )
}

