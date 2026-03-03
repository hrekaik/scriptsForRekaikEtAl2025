library(rtracklayer)
wd <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures/outputs/scRNAseq/integration"

dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

# Get the gene list from TOME github:
gene.list.name.file <- "mouse.v12.geneID.txt"
if (! file.exists(gene.list.name.file)){
  download.file(
    "https://raw.githubusercontent.com/ChengxiangQiu/tome_code/a41bed37cb87dc28fae0497cb7f4e6a7a91e33e1/help_code/mouse.v12.geneID.txt",
    gene.list.name.file
  )
}
gene.list.name <- read.delim(gene.list.name.file)

# # We use the an old version (86):
temp.file <- tempfile()
download.file(
  "https://ftp.ensembl.org/pub/release-86/gtf/mus_musculus/Mus_musculus.GRCm38.86.gtf.gz",
  temp.file
)
gtf <- readGFF(temp.file)
geneListName <- unique(gtf[, c("gene_id", "gene_name")])
colnames(geneListName) <- c("gene_ID", "gene_short_name_ensembl86")

# # We use the an old version (89):
temp.file <- tempfile()
download.file(
  "https://ftp.ensembl.org/pub/release-89/gtf/mus_musculus/Mus_musculus.GRCm38.89.gtf.gz",
  temp.file
)
gtf <- readGFF(temp.file)
geneListName89 <- unique(gtf[, c("gene_id", "gene_name")])
colnames(geneListName89) <- c("gene_ID", "gene_short_name_ensembl89")


# We use the genes we have in scRNAseq
gene.scrnaseq <- read.delim("../../../../toGEO/scRNAseq/Raw/GEX/HoxLess_CellPlex_1/genes.tsv", header = FALSE)
colnames(gene.scrnaseq) <- c("gene_ID", "gene_short_name_scrnaseq")
gene.scrnaseq$gene_short_name_scrnaseq <- make.unique(gene.scrnaseq$gene_short_name_scrnaseq)

all.genes <- merge(
  gene.list.name,
  merge(
    merge(
      geneListName, 
      geneListName89,
      by = "gene_ID", all = TRUE
      ),
    gene.scrnaseq,
    by = "gene_ID", all = TRUE
    ),
  by = "gene_ID", all = TRUE
  )

all.genes$nb <- apply(all.genes[, grep("gene_short_name", colnames(all.genes), value = TRUE)],
                      1,
                      function(v){length(unique(na.omit(v)))})
table(all.genes$nb)

#     1     2     3 
# 50879  4628  1032 



table(all.genes$nb, is.na(all.genes$gene_short_name_scrnaseq))
  #   FALSE  TRUE
  # 1 49539  1340
  # 2  4572    56
  # 3  1030     2



all.genes$name <- all.genes$gene_short_name_scrnaseq
all.genes$name[is.na(all.genes$gene_short_name_scrnaseq)] <- 
  all.genes$gene_short_name[is.na(all.genes$gene_short_name_scrnaseq)]
all.genes$name[is.na(all.genes$name)] <- 
  all.genes$gene_short_name_ensembl86[is.na(all.genes$name)]
all.genes$name[is.na(all.genes$name)] <- 
  all.genes$gene_short_name_ensembl89[is.na(all.genes$name)]

all.genes <- all.genes[order(all.genes$gene_short_name_scrnaseq), ]

all.genes$name <- make.unique(all.genes$name)

write.table(all.genes[, c("gene_ID", "name")], "gene_id_name_conversion.txt",
            row.names = FALSE, sep = "\t")
