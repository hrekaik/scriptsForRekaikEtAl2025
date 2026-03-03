library(GenomicRanges)
library(VennDiagram)
library(ggplot2)
library(grid)

wd <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures"
setwd(wd)
plot.dir <- "outputs/ChIPm"

wt_file <- "../toGEO/ChIPm/PBX1_macs2narrowPeak_wt.bed"
wt_bigwig <- "../toGEO/ChIPm/PBX1_wt.bigwig"
mutant_file <- "../toGEO/ChIPm/PBX1_macs2narrowPeak_BADC.bed"
mutant_bigwig <- "../toGEO/ChIPm/PBX1_BADC.bigwig"

pval_threshold <- 1e-5
fc_threshold <- 1

# Function to read narrowPeak files
read_narrowPeak <- function(file_path, sample_name) {
  peaks <- read.table(file_path, 
                      col.names = c("chr", "start", "end", "name", "score",
                                    "strand", "signalValue", "pValue", "qValue", "peak"))
  peaks <- peaks[peaks$pValue > -log10(pval_threshold), ]
  peaks$strand <- "*"
  gr <- GRanges(seqnames = peaks$chr,
                ranges = IRanges(start = peaks$start, end = peaks$end),
                strand = peaks$strand,
                score = peaks$score,
                signalValue = peaks$signalValue,
                pValue = peaks$pValue,
                qValue = peaks$qValue,
                peak = peaks$peak)
  mcols(gr)$sample <- sample_name
  return(gr)
}

# Read and process peaks
wt_peaks <- read_narrowPeak(wt_file, "wt")
mutant_peaks <- read_narrowPeak(mutant_file, "BADC")

# Merge overlapping peaks within each condition
merge_peaks <- function(gr) {
  merged <- GenomicRanges::reduce(gr, min.gapwidth = 0)
  hits <- findOverlaps(gr, merged)
  signal_vals <- tapply(gr$signalValue[queryHits(hits)], 
                        subjectHits(hits), mean)
  pvals <- tapply(gr$pValue[queryHits(hits)], 
                  subjectHits(hits), mean)
  mcols(merged)$signalValue <- signal_vals[as.character(1:length(merged))]
  mcols(merged)$pValue <- pvals[as.character(1:length(merged))]
  mcols(merged)$sample <- gr$sample[1]
  return(merged)
}
wt_merged <- merge_peaks(wt_peaks)
mutant_merged <- merge_peaks(mutant_peaks)

# Find overlapping and exclusive peaks
find_overlaps <- function(gr1, gr2, min_overlap = 1) {
  overlaps <- findOverlaps(gr1, gr2, minoverlap = min_overlap)
  gr1_overlap_idx <- unique(queryHits(overlaps))
  gr2_overlap_idx <- unique(subjectHits(overlaps))
  results <- list(
    gr1_only = gr1[-gr1_overlap_idx],
    gr2_only = gr2[-gr2_overlap_idx],
    common_gr1 = gr1[gr1_overlap_idx],
    common_gr2 = gr2[gr2_overlap_idx],
    common_count = length(gr1_overlap_idx)
  )
  return(results)
}
overlap_results <- find_overlaps(wt_merged, mutant_merged, min_overlap = 1)

# Create summary statistics
summary_stats <- data.frame(
  Category = c("WT only", "Mutant only", "Common"),
  Count = c(
    length(overlap_results$gr1_only),
    length(overlap_results$gr2_only),
    overlap_results$common_count
  )
)
print("Summary of peaks:")
print(summary_stats)

# Create Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(wt_merged),
  area2 = length(mutant_merged),
  cross.area = overlap_results$common_count,
  category = c("WT", "Mutant"),
  fill = c("#5f4b8b", "#e69a8d"),
  alpha = 0.7,
  cat.pos = c(0, 0),
  cat.dist = 0.05,
  scaled = TRUE
)
grid.newpage()
grid.draw(venn.plot)
ggsave(file.path(plot.dir, "PBX1_peaks_venn.pdf"), plot = venn.plot, width = 6, height = 6)




# Save common and exclusive peaks to BED files
export_granges <- function(gr, filename) {
  if (length(gr) == 0) {
    cat("Warning:", filename, "has 0 peaks. Creating empty file.\n")
    # Create an empty file with just a header comment
    writeLines("# Empty peak set - no peaks found", filename)
    return()
  }
  
  df <- data.frame(
    chr = as.character(seqnames(gr)),
    start = start(gr) - 1,  # BED is 0-based
    end = end(gr),
    name = paste0("peak_", 1:length(gr)),
    score = if ("score" %in% colnames(mcols(gr))) gr$score else rep(1000, length(gr)),
    strand = as.character(strand(gr))
  )
  df$strand[df$strand == "*"] <- "."
  write.table(df, filename, sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("Exported", length(gr), "peaks to", filename, "\n")
}

tryCatch({
  export_granges(overlap_results$gr1_only, file.path(plot.dir, "WT_exclusive_peaks.bed"))
}, error = function(e) {
  cat("Error exporting WT exclusive peaks:", e$message, "\n")
})

tryCatch({
  export_granges(overlap_results$gr2_only, file.path(plot.dir, "Mutant_exclusive_peaks.bed"))
}, error = function(e) {
  cat("Error exporting Mutant exclusive peaks:", e$message, "\n")
})

tryCatch({
  export_granges(overlap_results$common_gr1, file.path(plot.dir, "Common_peaks.bed"))
}, error = function(e) {
  cat("Error exporting Common peaks:", e$message, "\n")
})


export_granges(wt_merged, file.path(plot.dir, "WT_all_merged_peaks.bed"))
export_granges(mutant_merged, file.path(plot.dir, "Mutant_all_merged_peaks.bed"))

paste0(wd,"/",plot.dir,"/WT_exclusive_peaks.bed")

# Homer motifs
system2(
  "conda",
  args = c(
    "run", "-n", "pygenometracks",
    "findMotifsGenome.pl",
    paste0(wd,"/",plot.dir,"/Common_peaks.bed"),
    "mm10",
    paste0(wd,"/",plot.dir,"/Homer_peaks/Common_peaks"),
    "-size given -len 8,10,12")
)
system2(
  "conda",
  args = c(
    "run", "-n", "pygenometracks",
    "findMotifsGenome.pl",
    paste0(wd,"/",plot.dir,"/WT_exclusive_peaks.bed"),
    "mm10",
    paste0(wd,"/",plot.dir,"/Homer_peaks/WT_exclusive"),
    "-size given -len 8,10,12")
)
system2(
  "conda",
  args = c(
    "run", "-n", "pygenometracks",
    "findMotifsGenome.pl",
    paste0(wd,"/",plot.dir,"/Mutant_exclusive_peaks.bed"),
    "mm10",
    paste0(wd,"/",plot.dir,"/Homer_peaks/Mutant_exclusive"),
    "-size given -len 8,10,12")
)

### plotHeatmap
system2(
  "conda",
  args = c(
    "run", "-n", "pygenometracks",
    "computeMatrix", "reference-point",
    "--referencePoint", "center",
    "-b", "1000",
    "-a", "1000",
    "-R",
    paste0(wd,"/",plot.dir,"/WT_exclusive_peaks.bed"),
    paste0(wd,"/",plot.dir,"/Common_peaks.bed"),
    paste0(wd,"/",plot.dir,"/Mutant_exclusive_peaks.bed"),
    "-S",
    paste0(wd,"/",wt_bigwig),
    paste0(wd,"/",mutant_bigwig),
    "-o",
    paste0(wd,"/",plot.dir,"/matrix.gz"),
    "--missingDataAsZero"
  )
)

system2(
  "conda",
  args = c(
    "run",
    "-n",
    "pygenometracks",
    "plotHeatmap", 
    "-m",
    paste0(wd,"/",plot.dir,"/matrix.gz"),
    "-out",
    paste0(wd,"/",plot.dir,"/heatmap.pdf"),
    "--colorMap", "Spectral",
    "--zMin 0",
    "--zMax 10",
    "--whatToShow", "'heatmap and colorbar'",
    "--regionsLabel 'wt peaks' 'Common peaks' 'Hox-/- peaks'",
    "--samplesLabel wt Hox-/-",
    "--heatmapHeight 30"
    
  )
)

