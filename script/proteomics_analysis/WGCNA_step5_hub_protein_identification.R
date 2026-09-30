# ==============================================================================
# WGCNA Step 5: Hub Protein Identification
# Description: Calculate module membership (MM) and gene significance (GS),
#              identify hub proteins with high MM and high |GS|,
#              and visualize MM vs GS relationships
# ==============================================================================

library(WGCNA)
library(ggplot2)
library(dplyr)
library(svglite)

# ==============================================================================
# 1. Define paths and load data
# ==============================================================================
output_data <- "results/proteomics/WGCNA"
output_plot <- "figures/WGCNA"

load(file.path(output_data, "step4_module_trait_correlation.RData"))

# ==============================================================================
# 2. Define significant modules for analysis
# ==============================================================================
significant_modules_list <- c("turquoise", "black", "greenyellow", "blue", "brown")

# ==============================================================================
# 3. Calculate module membership (MM) and gene significance (GS)
# ==============================================================================
# Prepare trait data
datTraits <- data.frame(
  Familial = as.numeric(metadata_filtered$Group == "Familial")
)
rownames(datTraits) <- metadata_filtered$Sample_ID

# Calculate gene significance: correlation between each protein and trait
geneTraitCor <- cor(datExpr, datTraits$Familial, use = "p")
geneTraitPvalue <- corPvalueStudent(geneTraitCor, nrow(datExpr))

# Initialize results list
module_analysis_results <- list()

# Calculate MM and GS for each significant module
for(module_color in significant_modules_list) {
  
  # Extract proteins in this module
  module_proteins <- colnames(datExpr)[mergedColors == module_color]
  
  # Extract module eigengene
  module_column <- match(paste0("ME", module_color), colnames(MEs))
  module_eigengene <- MEs[, module_column]
  
  # Calculate MM: correlation between each protein and module eigengene
  MM <- cor(datExpr[, module_proteins], module_eigengene, use = "p")
  MMPvalue <- corPvalueStudent(MM, nrow(datExpr))
  
  # Extract GS for this module's proteins
  GS <- geneTraitCor[module_proteins, 1]
  GSPvalue <- geneTraitPvalue[module_proteins, 1]
  
  # Combine results
  module_df <- data.frame(
    Protein = module_proteins,
    Module = module_color,
    MM = as.vector(MM),
    MM_Pvalue = as.vector(MMPvalue),
    GS = GS,
    GS_Pvalue = GSPvalue,
    stringsAsFactors = FALSE
  )
  
  # Sort by MM (descending)
  module_df <- module_df[order(-abs(module_df$MM)), ]
  
  # Save to list
  module_analysis_results[[module_color]] <- module_df
}

# ==============================================================================
# 4. Identify hub proteins
# ==============================================================================
# Set thresholds for hub protein identification
MM_threshold <- 0.7
GS_threshold <- 0.3

hub_proteins_list <- list()

for(module_color in significant_modules_list) {
  
  module_df <- module_analysis_results[[module_color]]
  
  # Filter hub proteins: high MM and high |GS|
  hub_proteins <- module_df %>%
    filter(abs(MM) > MM_threshold & abs(GS) > GS_threshold) %>%
    arrange(desc(abs(MM)))
  
  hub_proteins_list[[module_color]] <- hub_proteins
}

# ==============================================================================
# 5. Save module protein lists and hub proteins
# ==============================================================================
# Save individual module results
for(module_color in significant_modules_list) {
  write.csv(module_analysis_results[[module_color]],
            file.path(output_data, sprintf("module_%s_proteins.csv", module_color)),
            row.names = FALSE)
  
  if(nrow(hub_proteins_list[[module_color]]) > 0) {
    write.csv(hub_proteins_list[[module_color]],
              file.path(output_data, sprintf("module_%s_hub_proteins.csv", module_color)),
              row.names = FALSE)
  }
}

# Save combined results
all_modules_df <- do.call(rbind, module_analysis_results)
write.csv(all_modules_df,
          file.path(output_data, "all_significant_modules_proteins.csv"),
          row.names = FALSE)

all_hubs_df <- do.call(rbind, hub_proteins_list)
if(nrow(all_hubs_df) > 0) {
  write.csv(all_hubs_df,
            file.path(output_data, "all_hub_proteins.csv"),
            row.names = FALSE)
}

# ==============================================================================
# 6. Visualize MM vs GS scatter plots
# ==============================================================================
# Combine all module data
plot_data <- do.call(rbind, module_analysis_results)

# Mark hub proteins
plot_data$IsHub <- abs(plot_data$MM) > MM_threshold & abs(plot_data$GS) > GS_threshold

# Create faceted scatter plot
p <- ggplot(plot_data, aes(x = MM, y = GS)) +
  geom_point(aes(color = IsHub), size = 1.5, alpha = 0.6) +
  geom_hline(yintercept = c(-GS_threshold, GS_threshold),
             linetype = "dashed", color = "gray50", size = 0.8) +
  geom_vline(xintercept = c(-MM_threshold, MM_threshold),
             linetype = "dashed", color = "gray50", size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", size = 0.8, alpha = 0.2) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "#E74C3C"),
                     labels = c("FALSE" = "Non-hub", "TRUE" = "Hub protein"),
                     name = "") +
  facet_wrap(~Module, scales = "free", ncol = 3) +
  labs(x = "Module Membership (MM)",
       y = "Gene Significance (GS) for Familial",
       title = "Module Membership vs Gene Significance") +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_plot, "WGCNA_step5_MM_vs_GS_scatter.png"),
       p, width = 12, height = 8, dpi = 600)

ggsave(file.path(output_plot, "WGCNA_step5_MM_vs_GS_scatter.svg"),
       p, width = 12, height = 8)

# ==============================================================================
# 7. Save intermediate data
# ==============================================================================
save(datExpr, metadata_filtered, mergedColors, MEs,
     module_analysis_results, hub_proteins_list,
     all_modules_df, all_hubs_df,
     MM_threshold, GS_threshold,
     file = file.path(output_data, "step5_module_analysis.RData"))

# ==============================================================================
# End of script
# ==============================================================================