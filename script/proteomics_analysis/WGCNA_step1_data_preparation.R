# ==============================================================================
# WGCNA Step 1: Data Preparation and Sample Clustering
# Description: Prepare expression matrix for WGCNA analysis, remove QC samples,
#              perform sample clustering to detect outliers
# ==============================================================================

library(WGCNA)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(svglite)

# WGCNA settings
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
set.seed(123)

# ==============================================================================
# 1. Define paths
# ==============================================================================
input_expr <- "results/proteomics/normalized_protein_matrix.csv"
input_metadata <- "data/proteomics/sample_metadata.csv"
output_plot <- "figures/WGCNA"
output_data <- "results/proteomics/WGCNA"

# Create output directories
dir.create(output_plot, recursive = TRUE, showWarnings = FALSE)
dir.create(output_data, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Load data
# ==============================================================================
expr_raw <- read.csv(input_expr, row.names = 1, check.names = FALSE)
metadata <- read.csv(input_metadata)

# ==============================================================================
# 3. Remove QC samples
# ==============================================================================
qc_cols <- grep("^QC", colnames(expr_raw), value = TRUE, ignore.case = TRUE)
total_int_cols <- grep("Total_Int", colnames(expr_raw), value = TRUE, ignore.case = TRUE)
exclude_cols <- c(qc_cols, total_int_cols)

expr_clean <- expr_raw[, !colnames(expr_raw) %in% exclude_cols]

# ==============================================================================
# 4. Align samples between expression matrix and metadata
# ==============================================================================
samples_in_expr <- colnames(expr_clean)
samples_in_metadata <- metadata$Sample_ID

# Keep common samples
common_samples <- intersect(samples_in_expr, samples_in_metadata)

metadata_filtered <- metadata[metadata$Sample_ID %in% common_samples, ]
metadata_filtered <- metadata_filtered[order(match(metadata_filtered$Sample_ID, common_samples)), ]
expr_clean <- expr_clean[, common_samples]

# Verify alignment
if(!all(colnames(expr_clean) == metadata_filtered$Sample_ID)) {
  stop("Error: Sample order mismatch between expression matrix and metadata")
}

# ==============================================================================
# 5. Prepare WGCNA input format (samples × proteins)
# ==============================================================================
datExpr <- as.data.frame(t(expr_clean))

# ==============================================================================
# 6. Sample clustering to detect outliers
# ==============================================================================
sampleTree <- hclust(dist(datExpr), method = "average")

# Plot sample dendrogram with group colors
traitColors <- numbers2colors(as.numeric(factor(metadata_filtered$Group)),
                              colors = c("blue", "red"))

# Save PNG
png(file.path(output_plot, "WGCNA_step1_sample_clustering.png"),
    width = 12, height = 8, units = "in", res = 600)
par(cex = 0.6, mar = c(4, 6, 2, 2))
plotDendroAndColors(sampleTree,
                    traitColors,
                    groupLabels = "Group",
                    main = "Sample Dendrogram and Group Colors")
dev.off()

# Save SVG
svglite(file.path(output_plot, "WGCNA_step1_sample_clustering.svg"),
        width = 12, height = 8)
par(cex = 0.6, mar = c(4, 6, 2, 2))
plotDendroAndColors(sampleTree,
                    traitColors,
                    groupLabels = "Group",
                    main = "Sample Dendrogram and Group Colors")
dev.off()

# ==============================================================================
# 7. Save intermediate data
# ==============================================================================
write.csv(datExpr,
          file.path(output_data, "datExpr_for_WGCNA.csv"),
          row.names = TRUE)

write.csv(metadata_filtered,
          file.path(output_data, "metadata_filtered.csv"),
          row.names = FALSE)

save(datExpr, metadata_filtered, sampleTree,
     file = file.path(output_data, "step1_data_preparation.RData"))

# ==============================================================================
# End of script
# ==============================================================================