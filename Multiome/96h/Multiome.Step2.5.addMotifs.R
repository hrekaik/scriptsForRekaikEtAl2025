library(devtools)
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("Signac")
# This package has the p6 version which we did not used
# safelyLoadAPackageInCRANorBioconductor("BSgenome.Mmusculus.UCSC.mm10")
safelyLoadAPackageInCRANorBioconductor("Rsamtools")
safelyLoadAPackageInCRANorBioconductor("BSgenome")
safelyLoadAPackageInCRANorBioconductor("TFBSTools") # To get motifs
safelyLoadAPackageInCRANorBioconductor("motifmatchr")
safelyLoadAPackageInCRANorBioconductor("JASPAR2022")
safelyLoadAPackageInCRANorBioconductor("chromVAR")

wd <- "/scratch/ldelisle/rstudio_test/"
setwd(wd)

genome.fasta <- "inputs_galaxy/mm10_UCSC.fa"
if (! file.exists(paste0(genome.fasta, ".fai"))) {
  indexFa(genome.fasta)
}
genome <- FaFile(genome.fasta)

combined.seurat <- readRDS("combined_Multiome_Velo.RDS")

# Add Motifs

# Get Motifs
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = TRUE)
)
names(pfm) <- paste0("JASPAR2022_", names(pfm))
# Get also from https://doi.org/10.1016/j.cell.2008.05.024
if (!file.exists("pwm_all_102107.txt")) {
  download.file(
    "https://hugheslab.ccbr.utoronto.ca/supplementary-data/homeodomains1/pwm_all_102107.txt",
    "pwm_all_102107.txt"
  )
}
lines.with.motifs <- readLines("pwm_all_102107.txt")
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
if (!file.exists("mm10_weight_matrices_v2")) {
  download.file(
    "https://swissregulon.unibas.ch/data/mm10_f5/mm10_weight_matrices_v2",
    "mm10_weight_matrices_v2"
  )
}
lines.with.motifs <- readLines("mm10_weight_matrices_v2")
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
  if (!file.exists(paste0("motif", i, ".motif"))) {
    download.file(
      paste0("http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults/motif", i, ".motif"),
      paste0("motif", i, ".motif")
    )
  }
}
pfm.from.paper <- list()
coef <- 1e6
for (i in seq_len(436)) {
  current.pfm <- pfm.from.homer(paste0("motif", i, ".motif"), coef = coef)
  pfm.from.paper[[ID(current.pfm)]] <- current.pfm
}
four.pfm <- c(three.pfm, pfm.from.paper)
table(sapply(lapply(four.pfm, Matrix), ncol))

saveRDS(four.pfm, "four.pfm.RDS")

DefaultAssay(combined.seurat) <- "ATAC"
motif <- AddMotifs(
  object = granges(combined.seurat),
  genome = genome,
  pfm = four.pfm
)

# See https://github.com/stuart-lab/signac/pull/1803

dim(motif)
# [1] 170070   2499
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

saveRDS(combined.seurat, "combined_Multiome_Velo_moreMotifs.RDS")
