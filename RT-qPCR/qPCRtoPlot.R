# Load necessary library
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
safelyLoadAPackageInCRANorBioconductor("tidyr")
safelyLoadAPackageInCRANorBioconductor("gridExtra")

filenames <- list.files(path = "RT-qPCR/data", pattern = "*.txt", full.names = TRUE)
plots <- list()
for (i in filenames) {
  gene <- gsub("*_wtvsBADC.txt", "", basename(i))
  data <- read.table(i, header = TRUE)
  data <- pivot_longer(data, cols = everything(), names_to = "Genotype", values_to = "deltaRps9")
  data$Genotype <- factor(data$Genotype, levels = c("wt", "BADC"))
  p <- ggplot(data, aes(x = Genotype, y = deltaRps9, color = Genotype)) +
    geom_jitter(width = 0.1, height = 0) +
    stat_compare_means(aes(group = Genotype),
      comparisons = list(c("wt", "BADC")),
      method = "t.test", paired = F,
      label = "p.signif"
    ) +
    labs(title = gene, y = "% Rps9") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.8, face = "italic"),
      axis.title.y = element_text(face = "italic")
    ) +
    scale_color_manual(values = c("#5f4b8b", "#e69a8d")) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)))

  plots[[gene]] <- p
}
dir.create("outputs/RT-qPCR/", showWarnings = FALSE, recursive = TRUE)
pdf(file = paste0("outputs/RT-qPCR/RT-qPCR.pdf"), width = 6, height = length(filenames) * 1.8)
set.seed(1)
grid.arrange(grobs = plots, ncol = 2)
dev.off()
