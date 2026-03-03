
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("Biostrings")
safelyLoadAPackageInCRANorBioconductor("BSgenome.Mmusculus.UCSC.mm10")
safelyLoadAPackageInCRANorBioconductor("GenomicRanges")
safelyLoadAPackageInCRANorBioconductor("ggplot2")


wd <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures"
setwd(wd)
plot.dir <- "outputs/Sangerseq"

fasta_file <- "Sangerseq/deletions_sequnce.fasta"
window_size <- 50

seqs <- readDNAStringSet(fasta_file)


clone_order <- names(seqs)

seq_df <- tibble(
  Sample = factor(clone_order, levels = clone_order),
  Sequence = as.character(seqs)
)


# Function to find center position
find_center_position <- function(seq) {
  insertion_match <- str_locate(seq, "-[ACGT]+-")
  if (!is.na(insertion_match[1])) {
    ins_start <- insertion_match[1] + 1
    ins_end   <- insertion_match[2] - 1
    return(floor((ins_start + ins_end) / 2))
  } else {
    return(str_locate(seq, "-")[1, 1])
  }
}


# Find center for each clone
seq_df <- seq_df %>%
  rowwise() %>%
  mutate(
    CenterPos = find_center_position(Sequence),
    Start = pmax(1, CenterPos - window_size),
    End   = pmin(nchar(Sequence), CenterPos + window_size),
    SubSeq = str_sub(Sequence, Start, End)
  ) %>%
  ungroup()


# Convert to long format
seq_long <- seq_df %>%
  dplyr::select(Sample, SubSeq) %>%
  mutate(SubSeq = strsplit(SubSeq, "")) %>%
  unnest(SubSeq) %>%
  group_by(Sample) %>%
  mutate(Position = row_number()) %>%
  ungroup()


# Plot
p <- ggplot(seq_long, aes(x = Position, y = Sample)) +
  geom_text(aes(label = SubSeq), family = "mono", size = 5) +
  scale_y_discrete(limits = rev(clone_order)) + 
  theme_minimal() +
  labs(
    x = "Position relative to event center",
    y = "Clone (FASTA order)",
    title = "Sequence context around deletion / insertion (±50 nt)"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )
ggsave(filename = file.path(plot.dir, "Plot_sequences.pdf"), plot = p, width = 14, height = 8)

