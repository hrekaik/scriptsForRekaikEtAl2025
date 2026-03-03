library(pheatmap)

gitHubDir <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures/RNAseq/RNAseq_ActA/ActAvsChrion"
output_dir <- file.path(gitHubDir,"../../../outputs/RNAseq/RNAseq_ActA/ActAvsChrion")

fpkms <- read.delim(file.path(gitHubDir,"../../../outputs/RNAseq/RNAseq_ActA/ActAvsChrion/AllCufflinks_Simplified.txt"))

samples.plans <- read.delim(samplesPlan<-file.path(gitHubDir,"samplesPlan_ActAvsChiron_wtonly_all.txt"))

gene_sets <- list(
  ActivinA_vs_Chiron = c(
    "Pou5f1", "Nanog", "Fgf5", "Otx2", "Cdx1", "Cdx2", "Cdx4",
    "Wnt3a", "Wnt5a", "Greb1", "Gdf11",
    "Wnt3", "Nodal", "Eomes", "T", "Fst",
    "Mixl1", "Tdgf1", "Lhx1", "Foxh1", "Pitx2",
    "Foxj1", "Noto", "Shh", "Nog", "Gsc",
    "Bmp7", "Eras", "Zic3", "Cer1", "Foxa2", "Sox17",
    "Gata4", "Gata6", "Dkk1", "Kdr", "Flt1",
    "Cdh11", "Mesp1"
  ),
  ActivinA_vs_Chiron_2 = c(
    "Fgf17",
    "Wnt8a", "Frzb",
    "Cdh1", "Epcam", "Krt8", "Krt18", "Cldn6",
    "Foxa1", "Cpm", "Foxd4",
    "Chrd", "Sox9",
    "Sox2", "Sox1", "Sox3", "Pax6", "Tbx6",
    "Dll1", "Mesp2", "Lfng", "Snai1",
    "Meox1", "Meox2", "Uncx", "Pax1", "Pax3", "Pax7",
    "Meis1", "Isl1", "Tbx1",
    "Gata1", "Tal1",
    "Lmx1a", "Fgf10", "Cd44", "Cd93", "Etv2",
    "Pecam1", "Aldh1a2", "Cyp26a1",
    "Meis2", "Lrp5", "Zic5", "Foxc2"
  ),
  Hox = c(
    ## Hoxa
    "Hoxa1", "Hoxa2", "Hoxa3", "Hoxa5", "Hoxa6", "Hoxa7",
    "Hoxa9", "Hoxa10", "Hoxa11", "Hoxa13",
    ## Hoxb
    "Hoxb1", "Hoxb2", "Hoxb3", "Hoxb4", "Hoxb5", "Hoxb6",
    "Hoxb7", "Hoxb8", "Hoxb9", "Hoxb13",
    ## Hoxc
    "Hoxc4", "Hoxc5", "Hoxc6", "Hoxc8", "Hoxc9",
    "Hoxc10", "Hoxc11", "Hoxc12", "Hoxc13",
    ## Hoxd
    "Hoxd1", "Hoxd3", "Hoxd4", "Hoxd8", "Hoxd9",
    "Hoxd10", "Hoxd11", "Hoxd12", "Hoxd13"
  )
)




rownames(samples.plans) <- samples.plans$sample
samples.plans$Line <- factor(samples.plans$Line,
                                  levels = c("ActA", "Chiron"))
samples.plans$Tissue <- factor(samples.plans$Tissue,
                               levels = c("96h", "120h"))
samples.plans <- samples.plans[
  order(samples.plans$Tissue,
        samples.plans$Line,
        samples.plans$Replicate),
]



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
    annotation_col = samples.plans[, c("Replicate", "Line", "Tissue")],
    annotation_colors = list(
      "Line" = c("ActA" = "#2bb673","Chiron" = "#ee4036"),
      "Replicate" = c("1" = "white", "2" = "grey70"),
      "Tissue" = c("96h" = "orange", "120h" = "#606c38")
    ),
    file = out_file,
    width = 6,
    height = 1+ (0.2 * nrow(mat)) # auto-scale height
  )
}

