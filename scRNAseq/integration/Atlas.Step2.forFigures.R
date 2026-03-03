library(usefulLDfunctions)

safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("reshape2")
safelyLoadAPackageInCRANorBioconductor("viridisLite")
safelyLoadAPackageInCRANorBioconductor("pheatmap")





wd <- "/Papers/DelHox/toGEO/scRNAseq"
dir.with.mouse.RDS <- "../../RDS_mouse_TOME"
integrated.rds <- file.path(dir.with.mouse.RDS, "seurat.atlas.combined.E6.5-9.5.gastruloids.integrated.study.rds")
fig.dir <- "../../Code_figures/outputs/scRNAseq/integration/"

setwd(wd)
dir.create(fig.dir, showWarnings = FALSE, recursive = TRUE)
mouse.HoxLess.gastruloid.integrated <- readRDS(integrated.rds)


original.Fates <- sort(na.omit(unique(mouse.HoxLess.gastruloid.integrated$Fate)))
mouse.HoxLess.gastruloid.integrated$Fate[is.na(mouse.HoxLess.gastruloid.integrated$Fate)] <- "NA"

short_name <- function(s) {
  paste(toupper(sapply(strsplit(gsub("\\(", "", s), " ")[[1]], substr, start = 1, stop = 1)), collapse = "")
}

cell_type_initials <- sapply(unique(mouse.HoxLess.gastruloid.integrated$cell_type), short_name)
cell_type_initials["Extraembryonic visceral endoderm"] <- "ExVE"
cell_type_initials["Splanchnic mesoderm"] <- "SpM"
cell_type_initials["Somatic mesoderm"] <- "SoM"
cell_type_initials["Epiblast"] <- "Ep"
cell_type_initials["Endothelium"] <- "En"
cell_type_initials["Retinal primordium"] <- "ReP"
cell_type_initials["Olfactory epithelium"] <- "OlE"
cell_type_initials["Otic epithelium"] <- "OtE"
cell_type_initials["Hindbrain"] <- "HB"
cell_type_initials["Neuromesodermal progenitors"] <- "NMPs"

mouse.HoxLess.gastruloid.integrated$short_cell_type <- unname(cell_type_initials[mouse.HoxLess.gastruloid.integrated$cell_type])


cat(paste0(sort(cell_type_initials), ": ", names(sort(cell_type_initials))), sep = ", ")
# A: Allantois, AER: Apical ectodermal ridge, AFP: Anterior floor plate,
# AM: Amniochorionic mesoderm, AMA: Amniochorionic mesoderm A,
# AMB: Amniochorionic mesoderm B, APS: Anterior primitive streak,
# BAE: Branchial arch epithelium, BP: Blood progenitors,
# CAOP: Chondrocyte and osteoblast progenitors, CLE: Caudal lateral epiblast,
# CN: Caudal neuroectoderm, D: Di/telencephalon, DE: Definitive endoderm,
# EE: Extraembryonic ectoderm, EM: Extraembryonic mesoderm, En: Endothelium,
# Ep: Epiblast, EVE: Embryonic visceral endoderm,
# ExVE: Extraembryonic visceral endoderm, F: Forebrain/midbrain,
# FE: Fusing epithelium, FHF: First heart field, G: Gut,
# GALE: Gut and lung epithelium, H: Hepatocytes, HB: Hindbrain,
# HP: Hematoendothelial progenitors, IM: Intermediate mesoderm,
# LMP: Limb mesenchyme progenitors, M: Mesencephalon/MHB, MM: Mixed mesoderm,
# MN: Motor neurons, MSC: Mesenchymal stromal cells, N: Notochord,
# NC: Neural crest, NCPG: Neural crest (PNS glia),
# NCPN: Neural crest (PNS neurons), NM: Nascent mesoderm,
# NMPs: Neuromesodermal progenitors, NPC: Neuron progenitor cells,
# OE: Olfactory epithelium, OPA: Osteoblast progenitors A,
# OPB: Osteoblast progenitors B, OtE: Otic epithelium, PA: Placodal area,
# PE: Parietal endoderm, PEC: Primitive erythroid cells,
# PFP: Posterior floor plate, PGC: Primordial germ cells,
# PK: Pre-epidermal keratinocytes, PMA: Paraxial mesoderm A,
# PMB: Paraxial mesoderm B, PMC: Paraxial mesoderm C,
# PSAAE: Primitive streak and adjacent ectoderm, RE: Renal epithelium,
# ReP: Retinal primordium, RN: Rostral neuroectoderm, RP: Roof plate,
# SC: Spinal cord, SCD: Spinal cord (dorsal), SCV: Spinal cord (ventral),
# SE: Surface ectoderm, SHF: Second heart field,
# SMP: Skeletal muscle progenitors, SoM: Somatic mesoderm,
# SpM: Splanchnic mesoderm

# Define colorcodes
my.colors <- 
  c(
    c("Unknown" = "#1f77b4", "Endothelial" = "#5254a3", "PGC" = "#e377c2",
      "Caudal Mesoderm" = "#ffbb78", "Caudal Epiblast" = "#e6550d",
      "Dermo/Sclero" = "#98df8a", "Somitic Mesoderm" = "#f7b6d2",
      "Anterior PSM" = "#ff9896", "Posterior PSM" = "#9467bd",
      "NMPs" = "#c5b0d5", "Neural Tube" = "#8c564b",
      "Mixed Mesoderm" = "#c49c94", "Nascent Mesoderm" = "#bcbd22",
      "Definitive Endoderm" = "#dbdb8d", "Primitive Streak" = "#1bca8d"),
    c("Caudal lateral epiblast" = "#758694", "Paraxial mesoderm B" = "#E0CCBE",
      "Paraxial mesoderm C" = "#AB886D",
      "Neuromesodermal progenitors" = "#c5b0d5",
      "Caudal neuroectoderm" = "#CDC2A5")
  )

default.color <- "grey"
not.defined <- na.omit(setdiff(mouse.HoxLess.gastruloid.integrated$cell_type, names(my.colors)))

my.colors.full <- c(rep(default.color, length(not.defined)), my.colors)
names(my.colors.full)[seq_along(not.defined)] <- not.defined

my.colors.full <- c(my.colors.full[names(cell_type_initials)], my.colors.full)
names(my.colors.full)[seq_along(cell_type_initials)] <- cell_type_initials


## figure for paper
tiff(file = paste0(fig.dir, "/GastruloidsvsMouse_E6.5-9.5.tiff"),
    width = 10 , height = 5, units = "in", res = 600, compression = "zip")
Idents(mouse.HoxLess.gastruloid.integrated) <- mouse.HoxLess.gastruloid.integrated$old.ident
p1 <- DimPlot(subset(mouse.HoxLess.gastruloid.integrated,subset = old.ident == "Mouse" | Genotype == "wt"),
              group.by  = "old.ident",
              reduction = "umap.cca.per.study",
              pt.size = 0.1,
              order = "Gastruloid",
              cols = c("#dcd5d5","#5f4b8b"),
              label.size = 4, raster=F) + theme(legend.position = "none")
p2 <- DimPlot(subset(mouse.HoxLess.gastruloid.integrated,subset = old.ident == "Mouse" | Genotype == "BADC"),
              group.by  = "old.ident",
              reduction = "umap.cca.per.study",
              pt.size = 0.1,
              order = "Gastruloid",
              cols = c("#dcd5d5","#e69a8d"),
              label.size = 4, raster=F) +
  theme(legend.position = "none")
print(ggarrange(p1, p2,
                ncol = 2, nrow = 1))
dev.off()
rm("p1", "p2")
gc()


xlow <- -4.5
xhigh <- 3
ylow <- -5.5
yhigh <- 0
x.poly <- c(3, 1, -2, -4.5, -2, 0)
y.poly <- c(-4.5, 0, -0.5, -3, -5, -5)

highlight <- function(g) {g + annotate("polygon", x = x.poly, y = y.poly,
                                       fill = NA, color = "black")}
zoom <- function(g) {highlight(g) + xlim(xlow, xhigh) + ylim(ylow, yhigh)}


tiff(file = paste0(fig.dir, "/GastruloidsvsMouse_E6.5-9.5_cluster.tiff"),
     width = 10 , height = 6, units = "in", res = 300)
Idents(mouse.HoxLess.gastruloid.integrated) <- mouse.HoxLess.gastruloid.integrated$old.ident
# temp <- subset(mouse.HoxLess.gastruloid.integrated, downsample = 2000)
gastruloid.df <- FetchData(subset(mouse.HoxLess.gastruloid.integrated, subset = old.ident == "Gastruloid"), vars = c("umapccaperstudy_1", "umapccaperstudy_2", "Genotype", "Fate"))
gastruloid.df <- gastruloid.df[sample(nrow(gastruloid.df)), ]
p1 <- DimPlot(subset(mouse.HoxLess.gastruloid.integrated,subset = old.ident == "Mouse"),
              group.by  = "Fate",
              reduction = "umap.cca.per.study",
              pt.size = 0.1,
              # order = c(rev(original.Fates), "NA"),
              cols = c("grey", "#ff9896","#e6550d","#ffbb78","#dbdb8d","#98df8a",
                       "#5254a3","#c49c94","#bcbd22","#8c564b","#c5b0d5",
                       "#e377c2","#9467bd","#1bca8d","#f7b6d2","#1f77b4"),
              label.size = 4, raster=F) +
  geom_point(data = subset(gastruloid.df, Genotype == "wt"), aes(x = umapccaperstudy_1, y = umapccaperstudy_2, color = Fate), size = 0.1) +
  # theme(legend.position = "none") +
  ggtitle("wt")
p2 <- DimPlot(subset(mouse.HoxLess.gastruloid.integrated,subset = old.ident == "Mouse"),
              group.by  = "Fate",
              reduction = "umap.cca.per.study",
              pt.size = 0.1,
              # order = c(rev(original.Fates), "NA"),
              cols = c("grey", "#ff9896","#e6550d","#ffbb78","#dbdb8d","#98df8a",
                       "#5254a3","#c49c94","#bcbd22","#8c564b","#c5b0d5",
                       "#e377c2","#9467bd","#1bca8d","#f7b6d2","#1f77b4"),
              label.size = 4, raster=F) +
  geom_point(data = subset(gastruloid.df, Genotype == "BADC"), aes(x = umapccaperstudy_1, y = umapccaperstudy_2, color = Fate), size = 0.1) +
  # theme(legend.position = "none") +
  ggtitle("BADC")
print(ggarrange(highlight(p1), highlight(p2),
                ncol = 2, nrow = 1, legend = "bottom", common.legend = TRUE))
dev.off()
rm("p1", "p2")
gc()

highlight(DimPlot(
  subset(mouse.HoxLess.gastruloid.integrated,subset = old.ident == "Mouse"),
  reduction = "umap.cca.per.study",
  group.by = "short_cell_type",
  label = TRUE
  ) +
  theme(legend.position = "none") +
    ggtitle("TOME E6.5 to E9.5"))
ggsave(file.path(fig.dir, "atlas_cell_types.pdf"), width = 15, height = 15)

highlight(DimPlot(
  subset(mouse.HoxLess.gastruloid.integrated,subset = old.ident == "Mouse"),
  reduction = "umap.cca.per.study",
  group.by = "day",
  cols = viridis(11),
  shuffle = TRUE
) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 3))) +
  theme(legend.position = "bottom") +
  ggtitle("TOME E6.5 to E9.5"))
ggsave(file.path(fig.dir, "atlas_days.pdf"), width = 15, height = 15)

sub.atlas <- subset(mouse.HoxLess.gastruloid.integrated, umapccaperstudy_1 < xhigh & umapccaperstudy_1 > xlow & umapccaperstudy_2 < yhigh & umapccaperstudy_2 > ylow)
table(sub.atlas$cell_type)
high.cell_type <- names(table(sub.atlas$cell_type))[table(sub.atlas$cell_type) > 500]
sub.atlas <- subset(sub.atlas, cell_type %in% high.cell_type)
sub.gastruloid.df <- subset(gastruloid.df, umapccaperstudy_1 < xhigh & umapccaperstudy_1 > xlow & umapccaperstudy_2 < yhigh & umapccaperstudy_2 > ylow)
# table(sub.gastruloid.df$Fate)
# high.cell_type <- names(table(sub.gastruloid.df$Fate))[table(sub.gastruloid.df$Fate) > 20]
# sub.gastruloid.df <- subset(sub.gastruloid.df, Fate %in% high.cell_type)

sub.sub.gastruloid.df <- subset(sub.gastruloid.df, Fate %in% c("Caudal Epiblast", "Caudal Mesoderm"))
# Get ellipse from my Fates
ellipses <- car::dataEllipse(sub.sub.gastruloid.df$umapccaperstudy_1, sub.sub.gastruloid.df$umapccaperstudy_2, levels=0.7, groups = factor(sub.sub.gastruloid.df$Fate))

add_ellipses <- function(g) {
  g +
    annotate(geom = "polygon", x = ellipses$`Caudal Epiblast`[, 1], y = ellipses$`Caudal Epiblast`[, 2], color = my.colors["Caudal Epiblast"], fill = NA, lty = 2, linewidth = 1) +
    annotate(geom = "polygon", x = ellipses$`Caudal Mesoderm`[, 1], y = ellipses$`Caudal Mesoderm`[, 2], color = my.colors["Caudal Mesoderm"], fill = NA, lty = 2, linewidth = 1)
    
}

g.sub.atlas.with.labels <- DimPlot(sub.atlas, reduction = "umap.cca.per.study", group.by = "short_cell_type", label = TRUE, pt.size = 0.5, raster = FALSE, shuffle = TRUE) +
  NoLegend() +
  ggtitle("TOME Atlas") +
  scale_color_manual(values = my.colors.full)
g.sub.atlas.with.day <- DimPlot(sub.atlas, reduction = "umap.cca.per.study", group.by = "day", label = TRUE, pt.size = 0.5, raster = FALSE, cols = viridis(11), shuffle = TRUE) +
  NoLegend() +
  ggtitle("TOME Atlas")
g.sub.wt <- ggplot(subset(sub.gastruloid.df, Genotype == "wt"), aes(x = umapccaperstudy_1, y = umapccaperstudy_2, color = Fate)) +
  geom_point(size = .5) +
  scale_color_manual(values = my.colors) +
  ggtitle("wt")  + 
  cowplot::theme_cowplot() +
  guides(color = guide_legend(override.aes = list(size = 2), nrow = 2)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom"
  )
g.sub.BADC <- ggplot(subset(sub.gastruloid.df, Genotype == "BADC"), aes(x = umapccaperstudy_1, y = umapccaperstudy_2, color = Fate)) +
  geom_point(size = .5) +
  scale_color_manual(values = my.colors) +
  ggtitle("BADC")  + 
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

ggarrange(
  ggarrange(add_ellipses(zoom(g.sub.atlas.with.day)), add_ellipses(zoom(g.sub.atlas.with.labels))),
  ggarrange(add_ellipses(zoom(g.sub.wt)), add_ellipses(zoom(g.sub.BADC)), common.legend = TRUE, legend = "bottom"),
  ncol = 1
)
ggsave(file.path(fig.dir, "Zoom_ellipses.pdf"), width = 10, height = 10)

# Get the legend for stages:
ggsave(file.path(fig.dir, "day_legend.pdf"),
       as_ggplot(get_legend(g.sub.atlas.with.day + theme(legend.position = "bottom"))),
       width = 5, height = 5)

# number/prop of cell per cluster and per day
num_cells <- mouse.HoxLess.gastruloid.integrated[[]] %>%
  filter(old.ident == "Mouse") %>%
  group_by(day, cell_type) %>%
  summarise(nb = n()) %>%
  mutate(prop = nb / sum(nb)) %>%
  filter(cell_type %in% c("Caudal lateral epiblast","Neuromesodermal progenitors","Paraxial mesoderm C"))

num_cells <- rbind(
  num_cells,
  data.frame(
    day = setdiff(unique(mouse.HoxLess.gastruloid.integrated$day, num_cells$day), c(NA, num_cells$day)),
    cell_type = "Caudal lateral epiblast",
    nb = 0,
    prop = 0
    )
  )

num_cells$cell_type <- factor(
  num_cells$cell_type,
  levels = c(
    "Caudal lateral epiblast",
    "Paraxial mesoderm C",
    "Neuromesodermal progenitors"
    )
)
ggplot(num_cells, aes(x = day, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(cell_type ~ ., labeller = labeller(cell_type = cell_type_initials)) +
  cowplot::theme_cowplot() +
  scale_y_continuous(
    name = "Percentage of cells per stage",
    labels = scales::percent
  ) +
  scale_fill_manual(values = my.colors.full)

ggsave(file.path(fig.dir, "Mouse_cluster_per_stage.pdf"), width = 9, height = 3.5)

## expression mouse vs Hox-/- gastruloids
mouse.HoxLess.gastruloid.integrated[["RNA"]] <- JoinLayers(mouse.HoxLess.gastruloid.integrated[["RNA"]])
gc()
library(scCustomize)
mouse.HoxLess.gastruloid.integrated$Genotype[is.na(mouse.HoxLess.gastruloid.integrated$Genotype)] <- "Mouse"
for (gene in c("Cdx1", "Grsf1", "Etv4", "Fgf3", "Fgf4", "Fgfbp3", "Cxcl12", "Igfbp5")) {
  tiff(
    file = file.path(fig.dir, paste0("expression_GastruloidsvsMouse_E6.5-9.5_", gene, ".tiff")),
    width = 3.5, height = 6, units = "in", res = 300
  )
  g <- FeaturePlot_scCustom(mouse.HoxLess.gastruloid.integrated,
    features = gene,
    split.by = "Genotype",
    alpha_exp = 0.8,
    pt.size = 0.3,
    num_columns = 1
  ) & NoAxes() & annotate("polygon",
    x = x.poly, y = y.poly,
    fill = NA, color = "black"
  )

  print(g)
  dev.off()
}

# Generate heatmaps
DefaultAssay(mouse.HoxLess.gastruloid.integrated) <- "RNA"
my.genes <- intersect(
  paste0("Hox", rep(letters[1:4], each = 13), rep(1:13, 4)),
  rownames(mouse.HoxLess.gastruloid.integrated)
)
# Select cells of interest:
mouse.HoxLess.gastruloid.integrated$global_stage <- ifelse(
  mouse.HoxLess.gastruloid.integrated$old.ident == "Mouse",
  mouse.HoxLess.gastruloid.integrated$day,
  mouse.HoxLess.gastruloid.integrated$Time
)
mouse.HoxLess.gastruloid.integrated$global_cell_type <- ifelse(
  mouse.HoxLess.gastruloid.integrated$old.ident == "Mouse",
  mouse.HoxLess.gastruloid.integrated$cell_type,
  mouse.HoxLess.gastruloid.integrated$Fate
)
mouse.HoxLess.gastruloid.integrated$global_suffix <- ifelse(
  mouse.HoxLess.gastruloid.integrated$old.ident == "Mouse",
  "",
  paste0("---", mouse.HoxLess.gastruloid.integrated$Genotype)
)
mouse.HoxLess.gastruloid.integrated$global_group <- paste0(
  mouse.HoxLess.gastruloid.integrated$global_cell_type, "---",
  mouse.HoxLess.gastruloid.integrated$global_stage,
  mouse.HoxLess.gastruloid.integrated$global_suffix
)
table(mouse.HoxLess.gastruloid.integrated$global_cell_type)
smaller.seurat <- subset(mouse.HoxLess.gastruloid.integrated, subset = global_cell_type %in% c("Neuromesodermal progenitors", "Caudal lateral epiblast", "NMPs", "Caudal Epiblast"))
with(smaller.seurat[[]], table(global_cell_type, global_stage))
#                              global_stage
# global_cell_type              120h  72h  96h E7.5 E7.75   E8 E8.25 E8.5a E9.5
#   Caudal Epiblast               15    8  479    0     0    0     0     0    0
#   Caudal lateral epiblast        0    0    0  380   576 1111   269     0    0
#   Neuromesodermal progenitors    0    0    0    0     0    0   985   963  265
#   NMPs                        1068    1   50    0     0    0     0     0    0
table.group <- table(smaller.seurat$global_group)
# Keep only groups with more than 50 cells
smaller.seurat <- subset(smaller.seurat, subset = global_group %in% names(table.group)[table.group > 50])

tmp <- PseudobulkExpression(
  object = smaller.seurat,
  group.by = "global_group",
  # features = my.genes, # I want to normalize to all genes so I must not set features here
  method = "aggregate",
  layer = "count",
  return.seurat = TRUE
)
matrix.average.norm <- tmp[["RNA"]]$data[my.genes, ]
annot.df <- unique(smaller.seurat[[]][, grep("global_", colnames(smaller.seurat[[]]))])
rownames(annot.df) <- annot.df$global_group
annot.df$global_group <- NULL

annot.df$Genotype <- NA
annot.df$Genotype[annot.df$global_suffix == "---wt"] <- "wt"
annot.df$Genotype[annot.df$global_suffix == "---BADC"] <- "BADC"
annot.df$Genotype <- factor(annot.df$Genotype, levels = c("wt", "BADC"))
annot.df$global_suffix <- NULL
annot.df$global_cell_type <- factor(
  annot.df$global_cell_type,
  levels = c("Caudal lateral epiblast", "Caudal Epiblast", "Neuromesodermal progenitors", "NMPs")
)
annot.df <- annot.df[order(annot.df$global_cell_type, annot.df$global_stage, annot.df$Genotype), ]
matrix.average.norm <- matrix.average.norm[, rownames(annot.df)]
library(pheatmap)
pheatmap(
  matrix.average.norm,
  annotation_col = annot.df,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_colors = list(
    global_cell_type = my.colors[as.character(unique(annot.df$global_cell_type))],
    Genotype = c("wt" = "#5f4b8b", "BADC" = "#e69a8d"),
    global_stage = c(
      "E7.5" = "#2A788EFF", "E7.75" = "#21908CFF", "E8" = "#22A884FF",
      "E8.25" = "#43BF71FF", "E8.5a" = "#7AD151FF", "E9.5" = "#BBDF27FF",
      "96h" = "#9CA986", "120h" = "#5F6F65"
    )
  ),
  gaps_col = c(4, 5, 8),
  filename = file.path(fig.dir, "Heatmap_Hox_NMP_CLE.pdf")
)

