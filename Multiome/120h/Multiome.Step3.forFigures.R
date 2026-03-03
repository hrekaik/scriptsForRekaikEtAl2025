if (!"usefulLDfunctions" %in% installed.packages()) {
  devtools::install_github("lldelisle/usefulLDfunctions")
}
if (!"presto" %in% installed.packages()) {
  devtools::install_github('immunogenomics/presto')
}

library(usefulLDfunctions)

safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("Signac")
safelyLoadAPackageInCRANorBioconductor("scCustomize")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("S4Vectors")
safelyLoadAPackageInCRANorBioconductor("GenomicRanges")
safelyLoadAPackageInCRANorBioconductor("patchwork")
safelyLoadAPackageInCRANorBioconductor("Matrix")
safelyLoadAPackageInCRANorBioconductor("stringr")
safelyLoadAPackageInCRANorBioconductor("tidyverse")
safelyLoadAPackageInCRANorBioconductor("viridis")



wd <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures"
setwd(wd)
plot.dir <- "outputs/Multiome/120h"


multiome.seurat <- readRDS("../toGEO/Multiome/120h/combined_Multiome_Velo_moreMotifs.RDS")
multiome.seurat$Genotype <- factor(x = multiome.seurat$Genotype, levels = c("wt", "BADC"))

# find clusters
multiome.seurat <- FindClusters(multiome.seurat, graph.name = "wknn")

### Change clusters name and order
Idents(multiome.seurat) <- multiome.seurat$seurat_clusters

DimPlot(multiome.seurat)

new.cluster.name <-  c(
  "0" = "Somitic_Mesoderm",
  "1" = "Somitic_Mesoderm",
  "2" = "NMPs",
  "3" = "Neural_tube",
  "4" = "Posterior_PSM",
  "5" = "Anterior_PSM",
  "6" = "Anterior_PSM",
  "7" = "Unknown",
  "8" = "PGCs",
  "9" = "Definitive_Endoderm"
)
my.colors <- c("#f7b6d2","#c5b0d5","#8c564b","#9467bd","#ff9896","#1f77b4","#e377c2","#dbdb8d")


multiome.seurat <- RenameIdents(multiome.seurat, new.cluster.name)
multiome.seurat[["Fate"]] <- Idents(multiome.seurat)


# build a UMAP visualization using RNA
multiome.seurat <- RunUMAP(
  object = multiome.seurat,
  reduction = "pca",
  assay = "SCT",
  dims = 1:30,
  verbose = TRUE
)

## UMAP by genotype
pdf(
  file = file.path(plot.dir, paste0("Multiome_bygeotype_RNAonly.pdf")),
  width = 5, height = 4
)
DimPlot(multiome.seurat,
        cols = c("#F5EDED", "#e69a8d"),
        label = T,
        pt.size = .1,
        alpha = 1,
        label.size = 4, group.by = "Genotype",
        order = c("BADC", "wt"),
        seed = 1
)
DimPlot(multiome.seurat,
        cols = c("#F5EDED", "#5f4b8b"),
        label = T,
        pt.size = .1,
        alpha = 1,
        label.size = 4, group.by = "Genotype",
        order = c("wt", "BADC"),
        seed = 1
)
dev.off()
p <- DimPlot(multiome.seurat, cols = my.colors)
ggsave(file.path(plot.dir, "Clustering_RNAonly.pdf"), p, width = 6.5, height = 4)

# build a joint UMAP visualization
multiome.seurat <- RunUMAP(
  object = multiome.seurat,
  nn.name = "weighted.nn",
  assay = "SCT",
  verbose = TRUE
)

## UMAP by genotype
pdf(
  file = file.path(plot.dir, paste0("Multiome_bygeotype.pdf")),
  width = 5, height = 4
)
DimPlot(multiome.seurat,
  cols = c("#F5EDED", "#e69a8d"),
  label = T,
  pt.size = .1,
  alpha = 1,
  label.size = 4, group.by = "Genotype",
  order = c("BADC", "wt"),
  seed = 1
)
DimPlot(multiome.seurat,
  cols = c("#F5EDED", "#5f4b8b"),
  label = T,
  pt.size = .1,
  alpha = 1,
  label.size = 4, group.by = "Genotype",
  order = c("wt", "BADC"),
  seed = 1
)
dev.off()




p <- DimPlot(multiome.seurat, cols = my.colors)
ggsave(file.path(plot.dir, "Clustering.pdf"), p, width = 6.5, height = 4)





## Marker dot plot
DefaultAssay(multiome.seurat) <- "SCT"
pdf(file = file.path(plot.dir, "dotplot.pdf"),
    width = 12 , height = 7.5)
DotPlot_scCustom(seurat_object = multiome.seurat, features = c("T", "Sox2","Nkx1-2", "Fgf8",
                                                               "Foxa2", "Sox17","Shh","Msgn1",
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


### NMP to Somitic mesoderm subset
caudal.tissue.seurat <- subset(multiome.seurat, cells = WhichCells(multiome.seurat,
                                                                   expression = Fate %in%  c("Somitic_Mesoderm", "Somitic_Mesoderm","NMPs",
                                                                                             "Posterior_PSM","Anterior_PSM","Anterior_PSM")))



### UMAP of selected genes
DefaultAssay(caudal.tissue.seurat) <- "SCT"
for (gene in c("Snai1","T", "Tbx6", "Foxc1", "Cyp26a1", "Hes7","Mesp2", "Ripply2", "Pax3","Cer1","Sox1","Sox2","Sox3",
               "Mesp1","Msgn1","Olig2","Twist2","Myf5","Bhlha15","Scx","Foxd3","Zic2","Lef1",
               "Pcdh8", "Pax1", "Pax7","Cdx1", "Grsf1", "Etv4", "Fgf3","Fgf4",
              "Hoxa1","Hoxa2", "Hoxa3", "Hoxb1","Hoxd1","Hoxb9",
               "Pbx1", "Meis1", "Pbx2", "Meis2")) {
  pdf(file = file.path(plot.dir,paste0("Multiome_Markers_", gene, ".pdf")),
      width = 3.5 , height = 6)
  g <- FeaturePlot_scCustom(caudal.tissue.seurat , features = gene,
                            split.by = "Genotype",
                            alpha_exp = 1,
                            pt.size = 0.8,
                            num_columns = 1) & NoAxes()
  print(g)
  dev.off()
}



#### pseudo time clustering heatmap for ATAC peaks

FeaturePlot_scCustom(caudal.tissue.seurat, features = "velocity_pseudotime")
ggsave(file.path(plot.dir, paste0("Pseudotime.pdf")), width = 6.5, height = 5)
caudal.tissue.seurat$pseudo.groupe <- cut(caudal.tissue.seurat$velocity_pseudotime, seq(0, 1, 0.05))
DimPlot(caudal.tissue.seurat, group.by = "pseudo.groupe",
        cols =  "Paired",
        pt.size = 1)
ggsave(file.path(plot.dir, paste0("pseudo.groupe.pdf")), width = 6.8, height = 5)

table(caudal.tissue.seurat$pseudo.groupe)
#(0,0.05] (0.05,0.1] (0.1,0.15] (0.15,0.2] (0.2,0.25] (0.25,0.3] (0.3,0.35] (0.35,0.4] (0.4,0.45] (0.45,0.5] (0.5,0.55] (0.55,0.6] (0.6,0.65] (0.65,0.7] (0.7,0.75] 
#0          0          0          0          0          0          0          4         34        559        457        153        315        326        342 
#(0.75,0.8] (0.8,0.85] (0.85,0.9] (0.9,0.95]   (0.95,1] 
#564       1445        813        746        402 

caudal.tissue.seurat$pseudo.groupe.genotype <- paste0(caudal.tissue.seurat$pseudo.groupe, "_", caudal.tissue.seurat$Genotype)
Idents(caudal.tissue.seurat) <- "pseudo.groupe.genotype"

table(caudal.tissue.seurat$pseudo.groupe.genotype)
#(0.35,0.4]_wt (0.4,0.45]_BADC   (0.4,0.45]_wt (0.45,0.5]_BADC   (0.45,0.5]_wt (0.5,0.55]_BADC   (0.5,0.55]_wt (0.55,0.6]_BADC   (0.55,0.6]_wt (0.6,0.65]_BADC 
#4               2              32             338             221             392              65             100              53             218 
#(0.6,0.65]_wt (0.65,0.7]_BADC   (0.65,0.7]_wt (0.7,0.75]_BADC   (0.7,0.75]_wt (0.75,0.8]_BADC   (0.75,0.8]_wt (0.8,0.85]_BADC   (0.8,0.85]_wt (0.85,0.9]_BADC 
#97             233              93             259              83             453             111            1080             365             212 
#(0.85,0.9]_wt (0.9,0.95]_BADC   (0.9,0.95]_wt   (0.95,1]_BADC     (0.95,1]_wt 
#601              13             733               2             400 



my.region.order <- c(
  "(0.45,0.5]_wt",
  "(0.5,0.55]_wt",
  "(0.55,0.6]_wt",
  "(0.6,0.65]_wt",
  "(0.65,0.7]_wt",
  "(0.7,0.75]_wt",
  "(0.75,0.8]_wt",
  "(0.8,0.85]_wt",
  "(0.85,0.9]_wt",
  "(0.45,0.5]_BADC",
  "(0.5,0.55]_BADC",
  "(0.55,0.6]_BADC",
  "(0.6,0.65]_BADC",
  "(0.65,0.7]_BADC",
  "(0.7,0.75]_BADC",
  "(0.75,0.8]_BADC",
  "(0.8,0.85]_BADC",
  "(0.85,0.9]_BADC"
)





if (!file.exists(file.path(plot.dir, "diffpeaks_allpseudogroupegenotype.csv"))) {
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


#####create dispersion plot
prepare_peaks <- function(clusters, Variable_peaks) {
  toppeaks <- Variable_peaks[Variable_peaks$cluster %in% clusters, ]
  toppeaks <- toppeaks[toppeaks$p_val_adj < 0.00001, ]
  toppeaks <- toppeaks %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
  
  toppeaks_filtered <- toppeaks %>%
    filter(stringr::str_count(gene, "-") < 3) %>%
    filter(!grepl("chrM", gene))
  length(toppeaks_filtered$gene)
  
  parsed_features <- do.call(rbind, strsplit(toppeaks_filtered$gene, "-"))
  colnames(parsed_features) <- c("chr", "start", "end")
  
  
  gr_df <- DataFrame(
    seqnames = parsed_features[, "chr"],
    start = as.numeric(parsed_features[, "start"]),
    end = as.numeric(parsed_features[, "end"]),
    strand = "*"
  )
  
  gr <- GRanges(
    seqnames = gr_df$seqnames,
    ranges = IRanges(start = gr_df$start, end = gr_df$end),
    strand = gr_df$strand
  )
  
  HoxY.coordinates <- GRanges(
    seqnames = c('chr6', "chr11", "chr15", "chr2", "chrY"),
    ranges = IRanges(start = c(52153796, 96191745, 102916767, 74667800,1),
                     end = c(52263152, 96372449, 103042756, 74767365,91744698)),
    strand = c('*', '*', '*', '*',"*")
  )
  HoxY.features.rows <- as.matrix(GenomicRanges::findOverlaps(gr,HoxY.coordinates))[,1]
  HoxY.features <- toppeaks_filtered$gene[HoxY.features.rows]
  toppeaks.minusHoxY <- setdiff(toppeaks_filtered$gene,HoxY.features)
}

calculate_pseudobulk <- function (seurat_object, genotype, peaks_list) {
  # pseudo bulk
  bulk <- PseudobulkExpression(
    object = subset(seurat_object, cells = WhichCells(seurat_object,
                                                      expression = Genotype == genotype)),
    assays = "ATAC",
    group.by = "pseudo.groupe.genotype",
    features = peaks_list,
    method = "average",
    layer = "counts"
  )
  bulk <- bulk[["ATAC"]][peaks_list, ]
}

calculate_interpolation <- function(mat_wt, mat_BADC) {
  ###### Interpolation
  
  # 1. Load  dgCMatrix and filter using most variable loci
  sparse_mat_wt <- mat_wt
  sparse_mat_BADC <- mat_BADC
  
  # 2. Extract time points from column names (using regex)
  time_intervals_wt <- str_match(colnames(sparse_mat_wt), "g\\((\\d+\\.\\d+),(\\d+\\.\\d+)\\]")
  time_points_wt <- as.numeric(time_intervals_wt[, 3])  # Use end of interval as time value
  
  time_intervals_BADC <- str_match(colnames(sparse_mat_BADC), "g\\((\\d+\\.\\d+),(\\d+\\.\\d+)\\]")
  time_points_BADC <- as.numeric(time_intervals_BADC[, 3])  # Use end of interval as time value
  
  # 3. Sort columns by time
  col_order_wt <- order(time_points_wt)
  sparse_mat_wt <- sparse_mat_wt[, col_order_wt]
  time_points_wt <- time_points_wt[col_order_wt]
  
  col_order_BADC <- order(time_points_BADC)
  sparse_mat_BADC <- sparse_mat_BADC[, col_order_BADC]
  time_points_BADC <- time_points_BADC[col_order_BADC]
  
  # 4. Define target time increments (e.g., 0.1)
  target_times <- seq(0.45, 0.9, by = 0.01)  # Adjust `from`, `to`, and `by` as needed
  
  # 5. Function to interpolate/extrapolate a single row
  interpolate_row_wt <- function(row_data) {
    approx(
      x = time_points_wt,
      y = row_data,
      xout = target_times,
      method = "linear",
      rule = 2  # Extrapolate outside the range
    )$y
  }
  
  interpolate_row_BADC <- function(row_data) {
    approx(
      x = time_points_BADC,
      y = row_data,
      xout = target_times,
      method = "linear",
      rule = 2  # Extrapolate outside the range
    )$y
  }
  
  # 6. Convert sparse matrix to dense (if feasible) and interpolate
  dense_mat_wt <- as.matrix(sparse_mat_wt)
  interp_mat_wt <- t(apply(dense_mat_wt, 1, interpolate_row_wt))
  
  dense_mat_BADC <- as.matrix(sparse_mat_BADC)
  interp_mat_BADC <- t(apply(dense_mat_BADC, 1, interpolate_row_BADC))
  
  # 7. Convert to data frame and set names
  interp_df_wt <- as.data.frame(interp_mat_wt)
  colnames(interp_df_wt) <- sprintf("t_%.2f", target_times)
  rownames(interp_df_wt) <- rownames(sparse_mat_wt)
  
  interp_df_BADC <- as.data.frame(interp_mat_BADC)
  colnames(interp_df_BADC) <- sprintf("t_%.2f", target_times)
  rownames(interp_df_BADC) <- rownames(sparse_mat_BADC)
  
  # 8. Combine wt and BADC data frames
  interp_df <- bind_rows(
    interp_df_wt %>% mutate(condition = "wt"),
    interp_df_BADC %>% mutate(condition = "BADC")
  )
  
}

peaks_dispersion <- function(df, group_id) {
  # Reshape to long format
  df_long <- df %>%
    mutate(region = rownames(df)) %>%
    pivot_longer(
      cols = starts_with("t_"),
      names_to = "time",
      names_prefix = "t_",
      names_transform = list(time = as.numeric),
      values_to = "value"
    ) %>%
    mutate(region = str_replace(region, "([^-]+-[^-]+-\\d+)[^0-9].*", "\\1"))
  
  # Find the timepoint of maximum value for each region and condition
  max_timepoints <- df_long %>%
    group_by(region, condition) %>%
    summarise(max_time = time[which.max(value)], .groups = 'drop')
  
  # Calculate differences and add group identifier
  final_df <- max_timepoints %>%
    pivot_wider(
      id_cols = region,
      names_from = condition,
      values_from = max_time
    ) %>%
    mutate(
      difference = wt - BADC,
      group_id = group_id  # Attach the group identifier
    )
  
  return(final_df)
}

all_dispersion_results <- list()

for (cluster in my.region.order) {
  time_intervals <- regmatches(cluster, 
                               regexec("\\(([0-9.]+),([0-9.]+)\\]_(wt|BADC)", cluster))[[1]]
  toppeaks.minusHoxY <- prepare_peaks(cluster,Variable_peaks_ATAC)
  
  bulk.pseudo.groupe.genotype.wt <- calculate_pseudobulk(caudal.tissue.seurat,"wt",toppeaks.minusHoxY)
  
  bulk.pseudo.groupe.genotype.BADC <- calculate_pseudobulk(caudal.tissue.seurat,"BADC",toppeaks.minusHoxY)
  
  interp_df <- calculate_interpolation (bulk.pseudo.groupe.genotype.wt, bulk.pseudo.groupe.genotype.BADC)
  
  # dispersion
  dispersion_result <- peaks_dispersion(interp_df,cluster)
  avg_diff <- mean(dispersion_result$difference, na.rm = TRUE)
  sd_value <- sd(dispersion_result$difference, na.rm = TRUE)
  all_dispersion_results[[cluster]] <- list(
    dispersion = dispersion_result,
    avg_diff = avg_diff,
    sd_value = sd_value
  )
}

#Combine dispersion results and average differences
combined_dispersion <- bind_rows(lapply(all_dispersion_results, function(x) x$dispersion))


combined_dispersion <- combined_dispersion %>%
  mutate(
    pseudotime = str_remove(group_id, "_.*$"),
    abs_diff   = abs(difference)
  )
#Create a dataframe of average differences
avg_diff_df <- combined_dispersion %>%
  group_by(pseudotime) %>%
  summarise(
    avg_abs_diff = mean(abs_diff, na.rm = TRUE),
    sd_value = sd(abs_diff, na.rm = TRUE),
    n = n(),                       
    .groups = "drop"
  )
avg_diff_df <- avg_diff_df %>%
  mutate(
    pt_start = as.numeric(str_extract(pseudotime, "(?<=\\()[^,]+")),
    pt_end   = as.numeric(str_extract(pseudotime, "(?<=,)[^\\]]+")),
    pt_mid   = (pt_start + pt_end) / 2
  )
p <- ggplot(
  data = avg_diff_df,
  aes(x = pt_mid, y = avg_abs_diff)
) +
  geom_line(color = "red", linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = avg_abs_diff - sd_value,
      ymax = avg_abs_diff + sd_value
    ),
    fill = "red",
    alpha = 0.2
  ) +
  labs(
    x = "Pseudotime",
    y = "Absolute Temporal shift between wt and Hox-/- peaks"
  ) +
  theme_classic()


ggsave(file.path(plot.dir, "Peaks_dispersion.pdf"), p, width = 5.2, height = 4, dpi = 300)



### Motif analysis over pseudotime
caudal.tissue.seurat <- subset(caudal.tissue.seurat, cells = WhichCells(caudal.tissue.seurat,
                                                                        expression = pseudo.groupe.genotype %in%  my.region.order))

motif_list <- c("HOXA1(Homeobox)/mES-Hoxa1-ChIP-Seq(SRP084292)/Homer",
                "Meis1(Homeobox)/MastCells-Meis1-ChIP-Seq(GSE48085)/Homer",
                "JASPAR2024-MA1113.1"
                
)
motif_list_names <- c("HOXA1(Homeobox)/mES-Hoxa1-ChIP-Seq(SRP084292)/Homer",
                      "Meis1(Homeobox)/MastCells-Meis1-ChIP-Seq(GSE48085)/Homer",
                      "JASPAR2024_MA1113.1"
)



chromvar_mat <- GetAssayData(
  caudal.tissue.seurat,
  assay = "chromvar",
  layer = "data"
)
for (n in seq_along(motif_list)) {
  motif_scores <- chromvar_mat[motif_list[n],]
  chromvar_df <- data.frame(
    cell = names(motif_scores),
    chromvar_score = as.numeric(motif_scores)
  )
  
  # Add metadata
  chromvar_df$pseudo.groupe <- caudal.tissue.seurat$pseudo.groupe
  chromvar_df$Genotype <- caudal.tissue.seurat$Genotype
  
  
  summary_data <- chromvar_df %>%
    group_by(pseudo.groupe, Genotype) %>%
    summarise(
      mean_signal = mean(chromvar_score, na.rm = TRUE),
      se_signal = sd(chromvar_score, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  summary_data <- summary_data %>%
    mutate(
      pseudotime_num = as.numeric(
        str_extract(as.character(pseudo.groupe), "(?<=,)[0-9.]+")
      )
    )
  
  p <- ggplot(
    summary_data,
    aes(
      x = pseudotime_num,
      y = mean_signal,
      color = Genotype,
      fill = Genotype,
      group = Genotype
    )
  ) +
    geom_ribbon(
      aes(
        ymin = mean_signal - se_signal,
        ymax = mean_signal + se_signal
      ),
      alpha = 0.25,
      color = NA
    ) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("wt" = "#5f4b8b", "BADC" = "#e69a8d")) +
    scale_fill_manual(values = c("wt" = "#5f4b8b", "BADC" = "#e69a8d")) +
    labs(
      title = paste(motif_list[n], "Motif Activity (chromVAR)"),
      x = "Pseudotime",
      y = "chromVAR Deviation Score",
      color = "Genotype",
      fill = "Genotype"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(
      limits = c(-4, 4),
      breaks = seq(-4, 4, by = 1)
    )
  
  
  ggsave(filename = paste0(plot.dir, "/",as.character(n), "_variation.pdf"), plot = p, width = 4, height = 4)
}





DefaultAssay(caudal.tissue.seurat) <- "ATAC"


for (n in seq_along(motif_list)) {
  p1 <- MotifPlot(caudal.tissue.seurat, motifs = motif_list_names[n])
  p2 <- FeaturePlot_scCustom(caudal.tissue.seurat, features = motif_list[n],
                             split.by = "Genotype",
                             colors_use = viridis(64, direction = -1)) &
    theme(legend.position = "none")
  p <- p1 +p2
  ggsave(filename = paste0(plot.dir, "/",as.character(n), "_logo.pdf"), plot = p, width = 12, height = 4)
}



### Finding overrepresented motifs with HOMER
toppeaks <- Variable_peaks_ATAC[Variable_peaks_ATAC$cluster %in% my.region.order, ]
toppeaks$cluster <- droplevels(toppeaks$cluster)
toppeaks <- toppeaks[toppeaks$p_val_adj < 0.00001, ]
toppeaks <- toppeaks %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
toppeaks <- toppeaks %>%
  separate(gene, into = c("seqnames", "start", "end"), sep = "-", convert = TRUE)

toppeaks.gr <- GenomicRanges::makeGRangesFromDataFrame(
  toppeaks,
  keep.extra.columns = TRUE
)

toppeaks.gr.split <- split(toppeaks.gr, toppeaks.gr$cluster)

for (cl in names(toppeaks.gr.split)) {
  
  gr <- GenomicRanges::sort(toppeaks.gr.split[[cl]])
  
  bed <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr) - 1,  # BED is 0-based
    end   = GenomicRanges::end(gr)
  )
  
  write.table(
    bed,
    file = file.path(plot.dir, paste0(cl, "_for_Homer.bed")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}


### Plot Homer motifs 
### Create a combined plot with both genotypes side by side

parse_motif_file <- function(file_path, folder_name) {
  # Extract pseudotime and genotype from folder name
  # Pattern: (start,end]_genotype_for_Homer_motifs
  folder_parts <- str_match(folder_name, "\\(([^,]+),([^]]+)\\]_([^_]+)_for_Homer_motifs")
  
  if (is.na(folder_parts[1,1])) {
    # Try alternative pattern without parentheses
    folder_parts <- str_match(folder_name, "([^,]+),([^]]+)_([^_]+)_for_Homer_motifs")
  }
  
  pseudotime <- paste0("(", folder_parts[1,2], ",", folder_parts[1,3], "]")
  genotype <- folder_parts[1,4]
  
  # Read file content
  content <- readLines(file_path)
  
  # Initialize results data frame
  results <- data.frame(
    Pseudotime = character(),
    Genotype = character(),
    Motif_name = character(),
    Score = numeric(),
    P_value = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each block starting with ">"
  current_block <- NULL
  
  for (line in content) {
    if (str_starts(line, ">")) {
      if (!is.null(current_block)) {
        parsed <- parse_motif_header(current_block, pseudotime, genotype)
        if (!is.null(parsed)) {
          results <- rbind(results, parsed)
        }
      }
      current_block <- line
    } else if (!is.null(current_block)) {
      if (str_detect(line, "^[0-9.\\t\\s-]+") && length(str_split(current_block, "\t")[[1]]) < 6) {
        current_block <- paste(current_block, line)
      }
    }
  }
  
  # Process the last block
  if (!is.null(current_block)) {
    parsed <- parse_motif_header(current_block, pseudotime, genotype)
    if (!is.null(parsed)) {
      results <- rbind(results, parsed)
    }
  }
  return(results)
}
parse_motif_header <- function(header_line, pseudotime, genotype) {
  parts <- str_split(header_line, "\t")[[1]]
  if (length(parts) < 6) {
    return(NULL)
  }
  
  # Extract BestGuess section (second part)
  best_guess_part <- parts[2]
  
  # Extract motif name and score
  motif_match <- str_match(best_guess_part, "BestGuess:([^/]+/[^/]+/[^/]+)\\(([^)]+)\\)")
  
  if (is.na(motif_match[1,1])) {
    motif_match <- str_match(best_guess_part, "BestGuess:([^(]+)\\(([^)]+)\\)")
  }
  
  if (is.na(motif_match[1,1])) {
    return(NULL)
  }
  
  motif_name <- motif_match[1,2]
  score <- as.numeric(motif_match[1,3])
  
  # Extract p-value (last part after P:)
  p_value_part <- parts[6]
  p_value_match <- str_match(p_value_part, "P:([^,]+)")
  
  if (is.na(p_value_match[1,1])) {
    if (str_detect(p_value_part, "P:")) {
      p_value <- str_extract(p_value_part, "(?<=P:)[^,]+")
    } else {
      p_value <- NA
    }
  } else {
    p_value <- p_value_match[1,2]
  }
  
  if (!is.na(p_value) && str_detect(p_value, "e-")) {
    p_value_str <- p_value
  } else if (!is.na(p_value)) {
    p_value_str <- p_value
  } else {
    p_value_str <- NA
  }
  
  return(data.frame(
    Pseudotime = pseudotime,
    Genotype = genotype,
    Motif_name = motif_name,
    Score = score,
    P_value = p_value_str,
    stringsAsFactors = FALSE
  ))
}

# Main function to process all directories
process_all_motif_directories <- function(root_dir = ".") {
  dirs <- list.dirs(root_dir, full.names = FALSE, recursive = FALSE)
  motif_dirs <- dirs[grepl("_for_Homer_motifs$", dirs)]
  
  if (length(motif_dirs) == 0) {
    stop("No directories matching the pattern '*_for_Homer_motifs' found.")
  }
  
  all_results <- data.frame(
    Pseudotime = character(),
    Genotype = character(),
    Motif_name = character(),
    Score = numeric(),
    P_value = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each directory
  for (dir in motif_dirs) {
    file_path <- file.path(root_dir, dir, "nonRedundant.motifs")
    
    if (file.exists(file_path)) {
      cat("Processing:", dir, "\n")
      dir_results <- parse_motif_file(file_path, dir)
      all_results <- rbind(all_results, dir_results)
    } else {
      cat("File not found:", file_path, "\n")
    }
  }
  
  return(all_results)
}

# Process all directories
motif_data <- process_all_motif_directories(plot.dir)

# View the results
print(head(motif_data))
print(dim(motif_data))

# Save to CSV
write.csv(motif_data, file.path(plot.dir, "motif_analysis_results.csv"), row.names = FALSE)


# Read the motif class mapping file
motif_class <- read.delim(file.path(plot.dir, "motif_class.txt"), sep = "\t", stringsAsFactors = FALSE)

# Clean column names
colnames(motif_class) <- trimws(colnames(motif_class))

# Remove duplicate motif entries (keep first occurrence)
motif_class <- motif_class %>%
  group_by(Motif_name) %>%
  dplyr::slice(1) %>%
  ungroup()

# Prepare the data with TF class mapping
motif_data_filtered <- motif_data %>%
  mutate(
    p_value_numeric = as.numeric(P_value),
    neg_log10_p = -log10(p_value_numeric)
  ) %>%
  filter(Score >= 0.7) %>%
  filter(neg_log10_p >= 12) %>%
  left_join(motif_class, by = "Motif_name") %>%
  filter(!is.na(TF) & TF != "UNKNOWN") %>%
  group_by(Pseudotime, Genotype, TF) %>%
  arrange(p_value_numeric, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(Pseudotime, Genotype) %>%
  arrange(p_value_numeric, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  group_by(Pseudotime, Genotype) %>%
  mutate(rank_in_group = row_number()) %>%
  ungroup()

write.csv(motif_data_filtered, file.path(plot.dir,"motif_analysis_filtered_results.csv"), row.names = FALSE)



# Define order of TF classes
tf_class_order <- c("NKX", "CDX", "TCF","TBX6", "MESP", "ZIC",  "FOX","HOX OCTA", "HOX DECA", "ETV", "TWIST", "NR")


create_combined_dot_plot <- function(data, tf_class_order = NULL, min_score = 0.7) {
  # Filter data for both genotypes
  plot_data <- data %>% 
    filter(Score >= min_score)
  
  if (nrow(plot_data) == 0) {
    message(paste("No motifs with score >=", min_score))
    return(NULL)
  }
  
  # Order pseudotime intervals 
  plot_data <- plot_data %>%
    mutate(
      pseudotime_start = as.numeric(gsub(".*\\(([0-9.]+),.*", "\\1", Pseudotime)),
      pseudotime_end = as.numeric(gsub(".*,([0-9.]+)\\].*", "\\1", Pseudotime)),
      pseudotime_mid = (pseudotime_start + pseudotime_end) / 2
    ) %>%
    arrange(pseudotime_mid) %>%
    mutate(Pseudotime = factor(Pseudotime, levels = unique(Pseudotime)))
  if (!is.null(tf_class_order)) {
    plot_data <- plot_data %>%
      filter(TF %in% tf_class_order) %>%
      mutate(TF = factor(TF, levels = rev(tf_class_order))) %>%
      mutate(Genotype = factor(Genotype, levels = c("wt","BADC")))
  } else {
    tf_class_order <- sort(unique(plot_data$TF))
    plot_data <- plot_data %>%
      mutate(TF = factor(TF, levels = rev(tf_class_order)))
  }
  
  # Create the combined plot with facets
  p <- ggplot(plot_data, aes(x = pseudotime_start, y = TF)) +
    
    geom_hline(
      yintercept = seq_along(unique(plot_data$TF)),
      color = "grey",
      linewidth = 0.3
    ) +
    
    geom_point(aes(size = Score, color = neg_log10_p), alpha = 0.8) +
    facet_grid(. ~ Genotype, scales = "free_x", space = "free") +
    
    scale_size_continuous(
      name = "Match score",
      range = c(2, 10),
      limits = c(min_score, 1),
      breaks = seq(min_score, 1, by = 0.05)
    ) +
    scale_color_gradientn(
      name = "-log10(P-value)",
      colors = plasma(64),
      limits = range(plot_data$neg_log10_p)
    ) +
    scale_x_continuous(
      breaks = seq(
        floor(min(plot_data$pseudotime_start)),
        ceiling(max(plot_data$pseudotime_start)),
        by = 0.05
      )
    ) +
    labs(
      title = paste("TF Class Enrichment (Score >=", min_score, ")"),
      x = "Pseudotime Interval",
      y = "TF Class"
    ) +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    guides(
      size = guide_legend(order = 1),
      color = guide_colorbar(order = 2)
    )
  
  
  return(p)
}

# Create combined plot with custom TF class order
combined_plot <- create_combined_dot_plot(motif_data_filtered, tf_class_order, min_score = 0.7)

if (!is.null(combined_plot)) {
  print(combined_plot)
  ggsave(file.path(plot.dir, "tf_class_combined_plot.pdf"), combined_plot, width = 10, height = 4.8, dpi = 300)
}