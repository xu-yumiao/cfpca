# ==============================================================================
# WGCNA Step 2: Soft Threshold Selection
# Description: Determine optimal soft-thresholding power for network construction
#              based on scale-free topology criterion
# ==============================================================================

library(WGCNA)
library(ggplot2)
library(svglite)
library(gridExtra)

# ==============================================================================
# 1. Define paths and load data from Step 1
# ==============================================================================
output_data <- "results/proteomics/WGCNA"
output_plot <- "figures/WGCNA"

load(file.path(output_data, "step1_data_preparation.RData"))

# ==============================================================================
# 2. Calculate network topology for different soft-thresholding powers
# ==============================================================================
powers <- c(1:10, seq(12, 20, by = 2))

sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         verbose = 5,
                         networkType = "signed")

# ==============================================================================
# 3. Extract and save fitting results
# ==============================================================================
sft_data <- sft$fitIndices

write.csv(sft_data,
          file.path(output_data, "soft_threshold_fitting.csv"),
          row.names = FALSE)

# ==============================================================================
# 4. Visualize scale-free topology fit
# ==============================================================================
plot_data <- data.frame(
  Power = sft_data$Power,
  SFT_R2 = -sign(sft_data$slope) * sft_data$SFT.R.sq,
  MeanK = sft_data$mean.k.
)

# Identify powers reaching R² thresholds
power_08 <- min(plot_data$Power[plot_data$SFT_R2 >= 0.8], na.rm = TRUE)
power_09 <- min(plot_data$Power[plot_data$SFT_R2 >= 0.9], na.rm = TRUE)
if(is.infinite(power_08)) power_08 <- NA
if(is.infinite(power_09)) power_09 <- NA

# Plot 1: Scale-free topology fit
p1 <- ggplot(plot_data, aes(x = Power, y = SFT_R2)) +
  geom_point(size = 3, color = "#2C3E50") +
  geom_line(size = 0.8, color = "#2C3E50") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#E74C3C", size = 0.8) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "#E67E22", size = 0.8) +
  geom_text(aes(label = Power), vjust = -0.5, size = 3) +
  annotate("text", x = max(powers) * 0.7, y = 0.82,
           label = "R² = 0.8", color = "#E74C3C", size = 3.5) +
  annotate("text", x = max(powers) * 0.7, y = 0.92,
           label = "R² = 0.9", color = "#E67E22", size = 3.5) +
  labs(x = "Soft Threshold (power)",
       y = "Scale Free Topology Model Fit, signed R²",
       title = "Scale Independence") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# Plot 2: Mean connectivity
p2 <- ggplot(plot_data, aes(x = Power, y = MeanK)) +
  geom_point(size = 3, color = "#2C3E50") +
  geom_line(size = 0.8, color = "#2C3E50") +
  geom_text(aes(label = Power), vjust = -0.5, size = 3) +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity",
       title = "Mean Connectivity") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# Combine plots
combined_plot <- grid.arrange(p1, p2, ncol = 2)

# Save plots
ggsave(file.path(output_plot, "WGCNA_step2_soft_threshold.png"),
       combined_plot,
       width = 10, height = 5, dpi = 600)

ggsave(file.path(output_plot, "WGCNA_step2_soft_threshold.svg"),
       combined_plot,
       width = 10, height = 5)

# ==============================================================================
# 5. Determine recommended soft threshold
# ==============================================================================
recommended_power <- sft$powerEstimate

# Fallback if automatic recommendation fails
if(is.na(recommended_power)) {
  if(!is.na(power_08)) {
    recommended_power <- power_08
  } else {
    recommended_power <- 6
  }
}

# ==============================================================================
# 6. Save results
# ==============================================================================
save(datExpr, metadata_filtered, sampleTree,
     sft, sft_data, recommended_power,
     file = file.path(output_data, "step2_soft_threshold.RData"))

# ==============================================================================
# End of script
# ==============================================================================

