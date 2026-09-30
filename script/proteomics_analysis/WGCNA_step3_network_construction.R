# ==============================================================================
# WGCNA Step 3: Network Construction and Module Detection
# Description: Build adjacency matrix and TOM, perform hierarchical clustering,
#              detect modules using dynamic tree cut, and merge similar modules
# ==============================================================================

library(WGCNA)
library(ggplot2)
library(svglite)

# ==============================================================================
# 1. Define paths and load data from Step 2
# ==============================================================================
output_data <- "results/proteomics/WGCNA"
output_plot <- "figures/WGCNA"

load(file.path(output_data, "step2_soft_threshold.RData"))

# ==============================================================================
# 2. Construct adjacency matrix and TOM
# ==============================================================================
# Calculate adjacency matrix using signed network
adjacency <- adjacency(datExpr,
                       power = recommended_power,
                       type = "signed")

# Convert to Topological Overlap Matrix (TOM)
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# ==============================================================================
# 3. Hierarchical clustering
# ==============================================================================
geneTree <- hclust(as.dist(dissTOM), method = "average")

# ==============================================================================
# 4. Module detection using dynamic tree cut
# ==============================================================================
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree,
                             distM = dissTOM,
                             deepSplit = 2,
                             pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

# Convert to color labels
dynamicColors <- labels2colors(dynamicMods)

# ==============================================================================
# 5. Calculate module eigengenes
# ==============================================================================
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

# Calculate dissimilarity between module eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

# ==============================================================================
# 6. Merge similar modules
# ==============================================================================
# Merge modules with correlation > 0.75 (dissimilarity < 0.25)
MEDissThres <- 0.25

# Plot module eigengene clustering tree
png(file.path(output_plot, "WGCNA_step3_module_merging_tree.png"),
    width = 10, height = 6, units = "in", res = 600)
plot(METree,
     main = "Clustering of Module Eigengenes",
     xlab = "",
     sub = "")
abline(h = MEDissThres, col = "red", lty = 2)
dev.off()

svglite(file.path(output_plot, "WGCNA_step3_module_merging_tree.svg"),
        width = 10, height = 6)
plot(METree,
     main = "Clustering of Module Eigengenes",
     xlab = "",
     sub = "")
abline(h = MEDissThres, col = "red", lty = 2)
dev.off()

# Merge modules
merge <- mergeCloseModules(datExpr,
                           dynamicColors,
                           cutHeight = MEDissThres,
                           verbose = 3)

mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# ==============================================================================
# 7. Visualize protein dendrogram and module colors
# ==============================================================================
png(file.path(output_plot, "WGCNA_step3_protein_dendrogram_modules.png"),
    width = 14, height = 8, units = "in", res = 600)
plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged Dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Protein Dendrogram and Module Colors")
dev.off()

svglite(file.path(output_plot, "WGCNA_step3_protein_dendrogram_modules.svg"),
        width = 14, height = 8)
plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged Dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Protein Dendrogram and Module Colors")
dev.off()

# ==============================================================================
# 8. Save module assignment
# ==============================================================================
moduleAssignment <- data.frame(
  Protein = colnames(datExpr),
  DynamicModule = dynamicColors,
  MergedModule = mergedColors,
  stringsAsFactors = FALSE
)

write.csv(moduleAssignment,
          file.path(output_data, "module_assignment.csv"),
          row.names = FALSE)

# ==============================================================================
# 9. Save intermediate data
# ==============================================================================
save(datExpr, metadata_filtered,
     adjacency, TOM, dissTOM, geneTree,
     dynamicMods, dynamicColors,
     mergedColors, mergedMEs,
     moduleAssignment,
     recommended_power,
     file = file.path(output_data, "step3_network_modules.RData"))

# ==============================================================================
# End of script
# ==============================================================================