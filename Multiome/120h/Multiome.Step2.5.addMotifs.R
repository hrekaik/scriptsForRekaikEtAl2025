library(devtools)
library(usefulLDfunctions)
remotes::install_github("stuart-lab/signac", ref = "develop")
safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("Signac")
# This package has the p6 version which we did not used
# safelyLoadAPackageInCRANorBioconductor("BSgenome.Mmusculus.UCSC.mm10")
safelyLoadAPackageInCRANorBioconductor("Rsamtools")
safelyLoadAPackageInCRANorBioconductor("BSgenome")

safelyLoadAPackageInCRANorBioconductor("TFBSTools") # To get motifs
safelyLoadAPackageInCRANorBioconductor("motifmatchr")
safelyLoadAPackageInCRANorBioconductor("JASPAR2024")
safelyLoadAPackageInCRANorBioconductor("chromVAR")

wd <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures"
setwd(wd)

genome.fasta <- "/Users/hocine.rekaik/Desktop/EXP/Gastruloids/scMultiome/21012025/data/mm10_UCSC.fa"
if (! file.exists(paste0(genome.fasta, ".fai"))) {
  indexFa(genome.fasta)
}
genome <- FaFile(genome.fasta)

combined.seurat <- readRDS("../toGEO/Multiome/120h/combined_Multiome_Velo.RDS")

# Add Motifs

# Get Motifs
# Get a list of motif position frequency matrices from the JASPAR database

# 1. Create the JASPAR2024 connection object
jaspar_obj <- JASPAR2024()

# 2. Extract the SQLite database connection from the object
# The 'db()' method accesses the database file path, which is then passed to RSQLite
jaspar_db_conn <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar_obj))

pfm <- getMatrixSet(
  x = jaspar_db_conn,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = TRUE)
)
names(pfm) <- paste0("JASPAR2024_", names(pfm))
# Get also from https://doi.org/10.1016/j.cell.2008.05.024
if (!file.exists("../toGEO/Multiome/120h/pwm_all_102107.txt")) {
  download.file(
    "https://hugheslab.ccbr.utoronto.ca/supplementary-data/homeodomains1/pwm_all_102107.txt",
    "data/pwm_all_102107.txt"
  )
}
lines.with.motifs <- readLines("../toGEO/Multiome/120h/pwm_all_102107.txt")
pfm.from.paper <- list()
matrix.values <- NULL
dimnames <- NULL
coef <- 1e6
for (i in seq_along(lines.with.motifs)) {
  if (i %% 7 == 1) {
    motif.name <- lines.with.motifs[i]
  }
  if (i %% 7 %in% c(3:6)) {
    letter <- strsplit(lines.with.motifs[i], ":\t")[[1]][1]
    dimnames <- c(dimnames, letter)
    values <- strsplit(gsub(paste0(letter, ":\t"), "", lines.with.motifs[i]), "\t")[[1]]
    matrix.values <- c(matrix.values, as.numeric(values))
  }
  if (i %% 7 == 6) {
    pfm.from.paper[[motif.name]] <- PFMatrix(ID = motif.name, name = motif.name,
                                             profileMatrix = matrix(coef * matrix.values,
                                                                    byrow = TRUE, nrow = 4,
                                                                    dimnames = list(dimnames)))
    matrix.values <- NULL
    dimnames <- NULL
  }
}
names(pfm.from.paper) <- paste0("Berger2008_", names(pfm.from.paper))
both.pfm <- c(pfm, pfm.from.paper)

# Get from ISMARA
if (!file.exists("../toGEO/Multiome/120h/mm10_weight_matrices_v2")) {
  download.file(
    "https://swissregulon.unibas.ch/data/mm10_f5/mm10_weight_matrices_v2",
    "data/mm10_weight_matrices_v2"
  )
}
lines.with.motifs <- readLines("../toGEO/Multiome/120h/mm10_weight_matrices_v2")
pfm.from.paper <- list()
matrix.values <- NULL
dimnames <- NULL
coef <- 1e2
for (my.line in lines.with.motifs) {
  if (substr(my.line, 1, 3) == "NA ") {
    motif.name <- gsub("^NA ", "", my.line)
  } else if (substr(my.line, 1, 2) == "P0") {
    dimnames <- strsplit(my.line, " ")[[1]]
    # Remove P0 and empty
    dimnames <- setdiff(dimnames, c("P0", ""))[1:4]
  } else if (my.line == "//") {
    pfm.from.paper[[motif.name]] <- PFMatrix(ID = motif.name, name = motif.name,
                                             profileMatrix = matrix(coef * matrix.values,
                                                                    byrow = FALSE, nrow = 4,
                                                                    dimnames = list(dimnames)))
    matrix.values <- NULL
  } else {
    values <- na.omit(as.numeric(strsplit(my.line, " ")[[1]]))
    pos <- as.numeric(values[1])
    stopifnot(length(matrix.values) == 4 * (pos - 1))
    matrix.values <- c(matrix.values, values[2:5])
  }
}
names(pfm.from.paper) <- paste0("ISMARA_", names(pfm.from.paper))
three.pfm <- c(both.pfm, pfm.from.paper)

# Get from HOMER
pfm.from.homer <- function(file.path, coef = 1e6) {
  lines.with.motifs <- readLines(file.path)
  motif.name <- strsplit(lines.with.motifs[1], "\t|,BestGuess")[[1]][2]
  matrix.values <- as.numeric(unlist(strsplit(lines.with.motifs[2:length(lines.with.motifs)], "\t")))
  return(
    PFMatrix(ID = motif.name, name = motif.name,
             profileMatrix = matrix(coef * matrix.values,
                                    byrow = FALSE, nrow = 4,
                                    dimnames = list(c("A", "C", "G", "T"))))
  )
  
}
for (i in seq_len(436)) {
  if (!file.exists(paste0("../toGEO/Multiome/120h/homer_motifs/motif", i, ".motif"))) {
    download.file(
      paste0("http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/motif", i, ".motif"),
      paste0("data/homer_motifs/motif", i, ".motif")
    )
  }
}
pfm.from.paper <- list()
coef <- 1e6
for (i in seq_len(436)) {
  current.pfm <- pfm.from.homer(paste0("../toGEO/Multiome/120h/homer_motifs/motif", i, ".motif"), coef = coef)
  pfm.from.paper[[ID(current.pfm)]] <- current.pfm
}
four.pfm <- c(three.pfm, pfm.from.paper)
table(sapply(lapply(four.pfm, Matrix), ncol))

saveRDS(four.pfm, "data/four.pfm.RDS")

DefaultAssay(combined.seurat) <- "ATAC"
# remove chr4-JH584295-random-1289-2046 fragment
combined.seurat[["ATAC"]] <- subset(combined.seurat[["ATAC"]],
                                    features = setdiff(rownames(combined.seurat[["ATAC"]]), "chr4-JH584295-random-1289-2046"))
#fragment.tokeep <- granges(combined.seurat)
#fragment.tokeep <- setdiff(fragment.tokeep,fragment.tokeep[124222])
motif <- AddMotifs(
  object = granges(combined.seurat),
  genome = genome,
  pfm = four.pfm
)

# See https://github.com/stuart-lab/signac/pull/1803

dim(motif)
# [1] 178793   2499
motif@data <- SeuratObject::CheckFeaturesNames(motif@data)
combined.seurat <- SetAssayData(
  object = combined.seurat,
  layer = "motifs",
  new.data = motif
)
dim(Motifs(combined.seurat))
# [1] 170070   2499

# Computing motif activities

# Set the number of CPU using
library(BiocParallel)
register(MulticoreParam(8)) # Use 8 cores

combined.seurat <- RunChromVAR(
  object = combined.seurat,
  genome = genome
)

saveRDS(combined.seurat, "../toGEO/Multiome/120h/combined_Multiome_Velo_moreMotifs.RDS")
