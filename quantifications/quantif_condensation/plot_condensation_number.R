# Load necessary library
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor('ggplot2')
safelyLoadAPackageInCRANorBioconductor('ggpubr')
safelyLoadAPackageInCRANorBioconductor('tidyr')

wd <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures"
setwd(wd)
plot.dir <- "outputs/quantifications/quantif_condensation"

data <- read.table("quantifications/quantif_condensation/quantification_both.txt", header = TRUE, fill = TRUE)
data <- pivot_longer(data, cols=everything(), names_to = "Genotype", values_to = "number")
data$Genotype <- factor(data$Genotype, levels = c("wt", "BADC"))

p <- ggplot(data, aes(x = Genotype, y = number)) +
  geom_boxplot(width = 0.5, aes(fill = Genotype)) +
  stat_compare_means(aes(group = Genotype),
    comparisons = list(c("wt", "BADC")),
    method = "t.test", paired = F,
    label = "p.signif"
  ) +
  labs(title = "Condensation number per gastruloid", y = "Condensation number per gastruloid") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.8),
    axis.title.y = element_text(face = "italic")
  ) +
  scale_fill_manual(values = c("#5F4B8BFF", "#E69A8DFF"))
pdf(file = file.path(plot.dir, "condensation_number.pdf"),width = 4, height = 5,useDingbats = FALSE)
print(p)
dev.off()
