# Detect large deletions using MACS2 peaks and bigWig coverage



safelyLoadAPackageInCRANorBioconductor("rtracklayer")
safelyLoadAPackageInCRANorBioconductor("GenomicRanges")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("zoo")
safelyLoadAPackageInCRANorBioconductor("GenomeInfoDb")


# Set working directory and input files

wd <- "/Papers/DelHox/toGEO/ATACseq"
plot.dir <- "../../Code_figures/outputs/ATACseq/"
setwd(wd)

bw_files <- list(
  wt1  = "mES_wt_ATAC_rep1.bigwig",
  wt2  = "mES_wt_ATAC_rep2.bigwig",
  mut1 = "mES_BADC_ATAC_rep1.bigwig",
  mut2 = "mES_BADC_ATAC_rep2.bigwig"
)

peak_files <- list(
  wt1  = "mES_wt_ATAC_rep1.narrowPeak.gz",
  wt2  = "mES_wt_ATAC_rep2.narrowPeak.gz",
  mut1 = "mES_BADC_ATAC_rep1.narrowPeak.gz",
  mut2 = "mES_BADC_ATAC_rep2.narrowPeak.gz"
)


# Import peaks and build anchor regions

peaks <- lapply(peak_files, import)
all_peaks <- GenomicRanges::reduce(unlist(GRangesList(peaks), use.names = FALSE))
anchor_regions <- resize(all_peaks, width = width(all_peaks) + 2000, fix = "center")


# Keep only canonical chromosomes and synchronize seqlevels

canonical_chr <- paste0("chr", c(1:22, "X", "Y"))
common_chr <- Reduce(
  intersect,
  list(
    canonical_chr,
    seqlevels(anchor_regions)
  )
)
common_chr <- canonical_chr[canonical_chr %in% common_chr]
anchor_regions <- keepSeqlevels(anchor_regions, common_chr, pruning.mode = "coarse")


# Import bigWigs
bw <- lapply(bw_files, function(f) import(f, as = "RleList"))
bw <- lapply(bw, function(x) x[common_chr])


# Extract mean ATAC signal per anchor region

get_signal <- function(bw_rle, regions) {
  binnedAverage(regions, bw_rle, "signal")$signal
}

signal <- list(
  wt1  = get_signal(bw$wt1,  anchor_regions),
  wt2  = get_signal(bw$wt2,  anchor_regions),
  mut1 = get_signal(bw$mut1, anchor_regions),
  mut2 = get_signal(bw$mut2, anchor_regions)
)

# Normalize signals and compute log2 ratios

norm <- function(x) x / median(x[x > 0])
signal_norm <- lapply(signal, norm)

sig_wt  <- rowMeans(cbind(signal_norm$wt1,  signal_norm$wt2))
sig_mut <- rowMeans(cbind(signal_norm$mut1, signal_norm$mut2))

df <- data.frame(
  chr = factor(as.character(seqnames(anchor_regions)), levels = common_chr),
  start = start(anchor_regions),
  end = end(anchor_regions),
  wt = sig_wt,
  mut = sig_mut,
  log2ratio = log2((sig_mut + 1e-6) / (sig_wt + 1e-6))
)


# Smooth log2 ratios along chromosomes

df <- df %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(smoothed = rollmedian(log2ratio, k = 5, fill = NA)) %>%
  ungroup()


# Prepare chromosome rectangles for plotting

chrom_info <- getChromInfoFromUCSC("mm10")
chrom_info <- chrom_info[chrom_info$chrom %in% common_chr, ]
chrom_info$chrom <- factor(chrom_info$chrom, levels = common_chr)

ymin <- floor(min(df$smoothed, na.rm = TRUE))
ymax <- ceiling(max(df$smoothed, na.rm = TRUE))

rect_df <- data.frame(
  chr = chrom_info$chrom,
  xmin = 0,
  xmax = chrom_info$size / 1e6,
  ymin = ymin - 1,
  ymax = ymax + 1
)


# Plot genome-wide signal with chromosome rectangles

p <- ggplot() +
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey90", alpha = 0.3) +
  geom_line(data = df, aes(x = start / 1e6, y = smoothed), linewidth = 0.6) +
  facet_grid(rows = vars(chr), scales = "free_x", switch = "y") +
  labs(
    x = "Genomic position (Mb)",
    y = "Smoothed log2(MUT / WT ATAC signal)",
    title = "Peak-anchored, bigWig-validated deletion scan"
  ) +
  theme_classic() +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    panel.spacing = unit(0.1, "lines")
  ) +
  scale_y_continuous(breaks = seq(ymin, ymax, by = 10))


ggsave(file.path(plot.dir, paste0("chromosomes_deletion_clean.pdf")), plot = p, width = 10, height = 12)
