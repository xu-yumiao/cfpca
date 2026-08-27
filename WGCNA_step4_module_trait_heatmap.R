# ==============================================================================
# WGCNA Step 4: Module-Trait Correlation Heatmap
# Description: Calculate and visualize correlations between module eigengenes
#              and clinical traits (Familial vs Sporadic)
# ==============================================================================

library(WGCNA)
library(ggplot2)
library(reshape2)
library(svglite)

# ==============================================================================
# 1. Define paths and load data
# ==============================================================================
output_data <- "results/proteomics/WGCNA"
output_plot <- "figures/WGCNA"

load(file.path(output_data, "step4_module_trait_correlation.RData"))

# ==============================================================================
# 2. Prepare trait data for both groups
# ==============================================================================
datTraits_both <- data.frame(
  Familial = as.numeric(metadata_filtered$Group == "Familial"),
  Sporadic = as.numeric(metadata_filtered$Group == "Sporadic")
)
rownames(datTraits_both) <- metadata_filtered$Sample_ID

# Calculate module-trait correlations
moduleTraitCor_both <- cor(MEs, datTraits_both, use = "p")
moduleTraitPvalue_both <- corPvalueStudent(moduleTraitCor_both, nrow(datExpr))

# ==============================================================================
# 3. Sort modules by correlation and significance
# ==============================================================================
module_data <- data.frame(
  Module = gsub("ME", "", rownames(moduleTraitCor_both)),
  Cor_Familial = moduleTraitCor_both[, "Familial"],
  Pval_Familial = moduleTraitPvalue_both[, "Familial"],
  Cor_Sporadic = moduleTraitCor_both[, "Sporadic"],
  Pval_Sporadic = moduleTraitPvalue_both[, "Sporadic"],
  stringsAsFactors = FALSE
)

# Classify modules by significance
sig_threshold <- 0.05
module_data$Group <- ifelse(module_data$Pval_Familial < sig_threshold & module_data$Cor_Familial > 0, "Positive",
                            ifelse(module_data$Pval_Familial < sig_threshold & module_data$Cor_Familial < 0, "Negative",
                                   "Non-significant"))

# Sort within groups
module_data$Order <- NA
module_data$Order[module_data$Group == "Positive"] <- rank(-module_data$Cor_Familial[module_data$Group == "Positive"])
module_data$Order[module_data$Group == "Non-significant"] <- rank(-abs(module_data$Cor_Familial[module_data$Group == "Non-significant"]))
module_data$Order[module_data$Group == "Negative"] <- rank(module_data$Cor_Familial[module_data$Group == "Negative"])

module_data$GroupOrder <- ifelse(module_data$Group == "Positive", 1,
                                 ifelse(module_data$Group == "Non-significant", 2, 3))

module_data <- module_data[order(module_data$GroupOrder, module_data$Order), ]
module_data$Module <- factor(module_data$Module, levels = module_data$Module)

# ==============================================================================
# 4. Prepare ordered matrices for heatmap
# ==============================================================================
cor_matrix <- moduleTraitCor_both
pval_matrix <- moduleTraitPvalue_both

cor_matrix_ordered <- cor_matrix[match(levels(module_data$Module), gsub("ME", "", rownames(cor_matrix))), ]
pval_matrix_ordered <- pval_matrix[match(levels(module_data$Module), gsub("ME", "", rownames(pval_matrix))), ]

rownames(cor_matrix_ordered) <- paste0("ME", levels(module_data$Module))
rownames(pval_matrix_ordered) <- paste0("ME", levels(module_data$Module))

# Prepare text labels with significance markers
textMatrix <- matrix("", nrow = nrow(cor_matrix_ordered), ncol = ncol(cor_matrix_ordered))
for(i in 1:nrow(cor_matrix_ordered)) {
  for(j in 1:ncol(cor_matrix_ordered)) {
    sig_mark <- ifelse(pval_matrix_ordered[i, j] < 0.001, "***",
                       ifelse(pval_matrix_ordered[i, j] < 0.01, "**",
                              ifelse(pval_matrix_ordered[i, j] < 0.05, "*", "")))
    textMatrix[i, j] <- paste0(sprintf("%.2f", cor_matrix_ordered[i, j]),
                               "\n", sig_mark)
  }
}

# ==============================================================================
# 5. Generate WGCNA-style heatmap
# ==============================================================================
png(file.path(output_plot, "WGCNA_step4_module_trait_heatmap.png"),
    width = 8, height = 10, units = "in", res = 600)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = cor_matrix_ordered,
               xLabels = colnames(cor_matrix_ordered),
               yLabels = rownames(cor_matrix_ordered),
               ySymbols = gsub("ME", "", rownames(cor_matrix_ordered)),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.2,
               cex.lab.x = 1.2,
               cex.lab.y = 1.0,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")
dev.off()

svglite(file.path(output_plot, "WGCNA_step4_module_trait_heatmap.svg"),
        width = 8, height = 10)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = cor_matrix_ordered,
               xLabels = colnames(cor_matrix_ordered),
               yLabels = rownames(cor_matrix_ordered),
               ySymbols = gsub("ME", "", rownames(cor_matrix_ordered)),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.2,
               cex.lab.x = 1.2,
               cex.lab.y = 1.0,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")
dev.off()

# ==============================================================================
# 6. Generate ggplot2-style heatmap (alternative visualization)
# ==============================================================================
heatmap_data <- melt(module_data[, c("Module", "Cor_Familial", "Cor_Sporadic")],
                     id.vars = "Module",
                     variable.name = "Trait",
                     value.name = "Correlation")

heatmap_data$Pvalue <- ifelse(heatmap_data$Trait == "Cor_Familial",
                              module_data$Pval_Familial[match(heatmap_data$Module, module_data$Module)],
                              module_data$Pval_Sporadic[match(heatmap_data$Module, module_data$Module)])

heatmap_data$Significance <- ifelse(heatmap_data$Pvalue < 0.001, "***",
                                    ifelse(heatmap_data$Pvalue < 0.01, "**",
                                           ifelse(heatmap_data$Pvalue < 0.05, "*", "")))

heatmap_data$Trait <- factor(heatmap_data$Trait,
                             levels = c("Cor_Familial", "Cor_Sporadic"),
                             labels = c("Familial", "Sporadic"))

p <- ggplot(heatmap_data, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Correlation)),
            size = 4.5, fontface = "bold", vjust = -0.3) +
  geom_text(aes(label = Significance),
            size = 5, vjust = 1.8, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Correlation") +
  labs(x = "", y = "",
       title = "Module-Trait Correlation") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

ggsave(file.path(output_plot, "WGCNA_step4_module_trait_heatmap_ggplot.png"),
       p, width = 7, height = 9, dpi = 600)

ggsave(file.path(output_plot, "WGCNA_step4_module_trait_heatmap_ggplot.svg"),
       p, width = 7, height = 9)

# ==============================================================================
# End of script
# ==============================================================================