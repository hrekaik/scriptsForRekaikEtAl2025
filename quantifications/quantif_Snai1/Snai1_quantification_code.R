safelyLoadAPackageInCRANorBioconductor("readxl")
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("stringr")
safelyLoadAPackageInCRANorBioconductor("purrr")
safelyLoadAPackageInCRANorBioconductor("ggplot2")

wd <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures"
setwd(wd)
plot.dir <- "outputs/quantifications/quantif_Snai1"


fluorescence_table <- read.csv("quantifications/quantif_Snai1/Fluorescence_Intensity_WT_BADC.csv")

#Normalize distance per sample
fluorescence_norm <- fluorescence_table %>%
  group_by(Sample) %>%
  mutate(
    Distance_norm = Distance_microns / max(Distance_microns, na.rm = TRUE),
    Condition = ifelse(str_detect(Sample, "wt"), "wt", "BADC")
  ) %>%
  ungroup()

# Min–max normalize intensity per sample
fluorescence_norm <- fluorescence_norm %>%  
  group_by(Sample) %>%
  mutate(
    Gray_Value_norm = (Gray_Value - min(Gray_Value, na.rm = TRUE)) /
      (max(Gray_Value, na.rm = TRUE) - min(Gray_Value, na.rm = TRUE))
  ) %>%
  ungroup()

# Bin normalized distance 
n_bins <- 100

fluorescence_binned <- fluorescence_norm %>%
  mutate(
    Distance_bin = cut(
      Distance_norm,
      breaks = seq(0, 1, length.out = n_bins + 1),
      include.lowest = TRUE,
      labels = FALSE
    )
  ) %>%
  group_by(Condition, Distance_bin) %>%
  summarise(
    Distance_norm = mean(Distance_norm, na.rm = TRUE),
    Mean_intensity = mean(Gray_Value_norm, na.rm = TRUE),
    SD_intensity = sd(Gray_Value_norm, na.rm = TRUE),
    .groups = "drop"
  )


# Summarize intensity per sample per region
regional_summary <- fluorescence_norm %>%
  mutate(
    Region = case_when(
      Distance_norm < 0.25  ~ "Anterior",
      Distance_norm >= 0.25 & Distance_norm < 0.6 ~ "DeterminationFront",
      Distance_norm >= 0.6 ~ "Posterior"
    )
  ) %>%
  group_by(Sample, Condition, Region) %>%
  summarise(
    Mean_intensity = mean(Gray_Value_norm, na.rm = TRUE),
    .groups = "drop"
  )
# Split data
anterior_data  <- regional_summary %>% filter(Region == "Anterior")
DeterminationFront_data <- regional_summary %>% filter(Region == "DeterminationFront")
posterior_data <- regional_summary %>% filter(Region == "Posterior")

# Wilcoxon tests
wilcox_anterior <- wilcox.test(
  Mean_intensity ~ Condition,
  data = anterior_data,
  exact = FALSE
)

wilcox_DeterminationFront <- wilcox.test(
  Mean_intensity ~ Condition,
  data = DeterminationFront_data,
  exact = FALSE
)

wilcox_posterior <- wilcox.test(
  Mean_intensity ~ Condition,
  data = posterior_data,
  exact = FALSE
)

wilcox_anterior
wilcox_DeterminationFront
wilcox_posterior




# Base plot with bars
regions_df <- data.frame(
  Region = c("Anterior", "DeterminationFront", "Posterior"),
  xmin = c(0, 0.25, 0.6),
  xmax = c(0.25, 0.6, 1),
  y = 1.05, 
  p_value = c(
    wilcox_anterior$p.value,
    wilcox_DeterminationFront$p.value,
    wilcox_posterior$p.value
  ),
  sig_label = c(
    ifelse(wilcox_anterior$p.value < 0.05, "*", "ns"),
    ifelse(wilcox_DeterminationFront$p.value < 0.05, "*", "ns"),
    ifelse(wilcox_posterior$p.value < 0.05, "*", "ns")
  )
)

p <- ggplot(fluorescence_binned,
       aes(x = Distance_norm, y = Mean_intensity,
           color = Condition, fill = Condition)) +
  
  geom_ribbon(
    aes(
      ymin = Mean_intensity - SD_intensity,
      ymax = Mean_intensity + SD_intensity
    ),
    alpha = 0.25,
    color = NA
  ) +
  scale_color_manual(values = c("wt" = "#5f4b8b", "BADC" = "#e69a8d")) +
  scale_fill_manual(values = c("wt" = "#5f4b8b", "BADC" = "#e69a8d")) +
  geom_line(size = 1.2) +
  geom_segment(
    data = regions_df,
    aes(x = xmin+0.01, xend = xmax, y = y, yend = y),
    inherit.aes = FALSE,
    color = "black",
    size = 1
  ) +
  geom_text(
    data = regions_df,
    aes(x = (xmin + xmax)/2, y = y + 0.05, label = sig_label),
    inherit.aes = FALSE,
    size = 5
  ) +
  
  labs(
    x = "Normalized Distance",
    y = "Fluorescence Intensity (Gray Value)",
    title = "Fluorescence Intensity Profile Along Normalized Distance"
  ) +
  
  scale_x_continuous(limits = c(0, 1)) +
  theme_classic()
ggsave(filename = file.path(plot.dir,"Snai1_wt_BADC_normalized_intensity.pdf"), plot = p, width = 12, height = 8)
