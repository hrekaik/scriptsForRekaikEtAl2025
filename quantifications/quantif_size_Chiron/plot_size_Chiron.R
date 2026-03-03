
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor('ggplot2')
safelyLoadAPackageInCRANorBioconductor('ggpubr')
safelyLoadAPackageInCRANorBioconductor('tidyr')
safelyLoadAPackageInCRANorBioconductor('dplyr')

set.seed(123) 

wd <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures"
setwd(wd)
plot.dir <- "outputs/quantifications/quantif_size_Chiron"

# Read data
data <- read.table("quantifications/quantif_size_Chiron/quantif_size.txt", header = TRUE, fill = TRUE)

# Reshape + extract genotype and replicate
data_long <- data %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Len"
  ) %>%
  mutate(
    Genotype = case_when(
      grepl("^wt", Sample)   ~ "wt",
      grepl("^BADC", Sample) ~ "BADC"
    ),
    Rep = case_when(
      grepl("rep1", Sample) ~ "rep1",
      grepl("rep2", Sample) ~ "rep2",
      grepl("rep3", Sample) ~ "rep3",
      grepl("rep4", Sample) ~ "rep4",
      grepl("rep5", Sample) ~ "rep5"
    )
  ) %>%
  filter(!is.na(Len))

# Subset 20 values per Genotype × Rep and compute mean
means_rep <- data_long %>%
  group_by(Genotype, Rep) %>%
  slice_sample(n = 20) %>%        
  summarise(
    Mean_Len = mean(Len),
    .groups = "drop"
  )

# Set factor order
means_rep$Genotype <- factor(means_rep$Genotype, levels = c("wt", "BADC"))

# percentage decrease
summary_df <- means_rep %>%
  group_by(Genotype) %>%
  summarise(mean_len = mean(Mean_Len, na.rm = TRUE))

wt_mean <- summary_df$mean_len[summary_df$Genotype == "wt"]
mut_mean <- summary_df$mean_len[summary_df$Genotype == "BADC"]

percent_decrease <- (wt_mean - mut_mean) / wt_mean * 100

percent_decrease

p <- ggplot(means_rep, aes(x = Genotype, y = Mean_Len)) +
  geom_point(
    aes(shape = Rep, color = Genotype),
    size = 8,
    position = position_jitter(width = 0.1)
  ) +
  ylim(0, 1000) +
  stat_compare_means(
    aes(group = Genotype),
    method = "t.test",
    paired = FALSE,
    label = "p.signif"
  ) +
  scale_color_manual(values = c("#5F4B8BFF", "#E69A8DFF")) +
  labs(
    title = "Length (replicate means)",
    y = "Mean length (µm)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.8, face = "italic"),
    axis.title.y = element_text(face = "italic")
  )



pdf(file = file.path(plot.dir, "size_Chiron.pdf"), width = 4, height = 5, useDingbats = FALSE)
p
dev.off()
