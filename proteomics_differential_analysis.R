# ==============================================================================
# Differential Protein Expression Analysis using Limma
# Description: Identify differentially expressed proteins (DEPs) between 
#              Familial and Sporadic groups using limma with empirical Bayes
# Criteria: |FC| > 1.5, adjusted p-value < 0.05
# ==============================================================================

library(limma)
library(dplyr)

# ==============================================================================
# 1. Define paths
# ==============================================================================
matrix_path <- "results/proteomics/normalized_protein_matrix.csv"
meta_path   <- "data/proteomics/sample_metadata.csv"
output_dir  <- "results/proteomics"

# ==============================================================================
# 2. Load normalized expression matrix and metadata
# ==============================================================================
expr_matrix <- read.csv(matrix_path, row.names = 1, check.names = FALSE)
metadata <- read.csv(meta_path, stringsAsFactors = FALSE)

# ==============================================================================
# 3. Remove QC samples
# ==============================================================================
non_qc_idx <- !grepl("QC", colnames(expr_matrix), ignore.case = TRUE)
expr_matrix_clean <- expr_matrix[, non_qc_idx]

# Align metadata with expression matrix
metadata_clean <- metadata[match(colnames(expr_matrix_clean), metadata$Sample_ID), ]

# ==============================================================================
# 4. Build design matrix for Familial vs Sporadic comparison
# ==============================================================================
valid_idx <- !is.na(metadata_clean$Group)
expr_matrix_clean <- expr_matrix_clean[, valid_idx]
group_factor <- factor(metadata_clean$Group[valid_idx])

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)

# Define contrast
contrast.matrix <- makeContrasts(
  contrasts = "Familial - Sporadic", 
  levels = design
)

# ==============================================================================
# 5. Perform limma differential analysis
# ==============================================================================
fit <- lmFit(expr_matrix_clean, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract results for all proteins
res_all <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
res_all$Protein <- rownames(res_all)

# ==============================================================================
# 6. Classify DEPs based on significance criteria
# ==============================================================================
res_all <- res_all %>%
  mutate(
    Significance = case_when(
      adj.P.Val < 0.05 & logFC > log2(1.5) ~ "Up",
      adj.P.Val < 0.05 & logFC < -log2(1.5) ~ "Down",
      TRUE ~ "Not Sig"
    )
  )

# Save results
write.csv(res_all, 
          file.path(output_dir, "differential_proteins_results.csv"), 
          row.names = FALSE, quote = FALSE)

# ==============================================================================
# End of script
# ==============================================================================