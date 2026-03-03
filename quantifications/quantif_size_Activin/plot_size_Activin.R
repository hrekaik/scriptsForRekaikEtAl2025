library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor('ggplot2')
safelyLoadAPackageInCRANorBioconductor('ggpubr')
safelyLoadAPackageInCRANorBioconductor('tidyr')

wd <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures"
setwd(wd)
plot.dir <- "outputs/quantifications/quantif_size_Activin"

data <- read.table("quantifications/quantif_size_Activin/quantification_activinA.txt",header = TRUE, fill = TRUE)
data <- pivot_longer(data, cols=everything(), names_to = "Genotype", values_to = "Len")
data$Genotype <- factor(data$Genotype, levels = c("wt","BADC"))
p <- ggplot(data, aes(x = Genotype, y = Len)) +
  geom_violin(aes(fill = Genotype)) +
  geom_boxplot(width = 0.2, alpha=0.7) +
  ylim(0,400) +
  stat_compare_means(aes(group = Genotype),method = "t.test",paired = F,
                     label = "p.signif") +
  labs(title = "Length", y="Length (µm)") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.8,face="italic"),
        axis.title.y=element_text(face="italic")) +
  scale_fill_manual(values=c("#5F4B8BFF","#E69A8DFF"))
pdf(file = file.path(plot.dir, "size_Activin.pdf"),width = 4, height = 5,useDingbats = FALSE)
p
dev.off()

