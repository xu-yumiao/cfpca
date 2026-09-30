# ==============================================================================
# Proteomics Data Preprocessing and Quantile Normalization
# Description: Process raw DIA proteomics data including filtering, 
#              log2 transformation, KNN imputation, quantile normalization, 
#              and PCA visualization
# ==============================================================================

# Load required libraries
library(readxl)
library(dplyr)
library(limma)
library(impute)
library(preprocessCore)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(svglite)

# ==============================================================================
# 1. Define paths (adjust according to your directory structure)
# ==============================================================================
input_dir <- "data/proteomics/raw"
mapping_path <- "data/proteomics/sample_mapping.xlsx"
output_dir <- "results/proteomics"
output_plot <- "figures"

# Create output directories if they don't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_plot, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Load raw protein data
# ==============================================================================
raw_data <- read_excel(file.path(input_dir, "List1_Proteins.xlsx"))

# Filter valid gene names and remove duplicates based on total intensity
data_clean <- raw_data %>% 
  filter(!is.na(`Gene Name`) & `Gene Name` != "") %>% 
  mutate(Total_Int = rowSums(select(., 10:ncol(.)), na.rm = TRUE)) %>% 
  arrange(desc(Total_Int)) %>% 
  distinct(`Gene Name`, .keep_all = TRUE)

# Extract expression matrix
expr_matrix_raw <- as.data.frame(lapply(data_clean %>% select(10:ncol(.)), as.numeric)) 
rownames(expr_matrix_raw) <- data_clean$`Gene Name`

# ==============================================================================
# 3. Sample ID mapping
# ==============================================================================
mapping_df <- read_excel(mapping_path)
meta_info <- mapping_df %>% 
  select(Group = 2, Sample_ID = 3, Report_Name = 4) %>% 
  filter(!is.na(Report_Name)) %>%
  mutate(Group = trimws(Group), 
         Patient_ID = sub("-[0-9]+$", "", Sample_ID))

# Replace column names with actual sample IDs
match_idx <- match(colnames(expr_matrix_raw), meta_info$Report_Name)
new_colnames <- meta_info$Sample_ID[match_idx]
new_colnames[is.na(new_colnames)] <- colnames(expr_matrix_raw)[is.na(new_colnames)]
colnames(expr_matrix_raw) <- new_colnames

# ==============================================================================
# 4. Log2 transformation and missing value filtering
# ==============================================================================
expr_matrix_raw[expr_matrix_raw == 0] <- NA
expr_log2 <- log2(expr_matrix_raw)

# Remove proteins with >50% missing values
keep_idx <- rowSums(is.na(expr_log2)) < (ncol(expr_log2) * 0.5)
expr_filtered <- expr_log2[keep_idx, ]

# ==============================================================================
# 5. KNN imputation and quantile normalization
# ==============================================================================
# KNN imputation for missing values
capture.output({ 
  expr_imputed <- impute.knn(as.matrix(expr_filtered), k = 10)$data 
})

# Quantile normalization
expr_norm_q <- normalize.quantiles(as.matrix(expr_imputed))
rownames(expr_norm_q) <- rownames(expr_imputed)
colnames(expr_norm_q) <- colnames(expr_imputed)
expr_final_q <- as.data.frame(expr_norm_q)

# Save normalized matrix
write.csv(expr_final_q, 
          file.path(output_dir, "normalized_protein_matrix.csv"), 
          quote = FALSE)

# ==============================================================================
# 6. PCA analysis and visualization
# ==============================================================================
# Prepare metadata for all samples
all_samples <- colnames(expr_final_q)
full_meta <- data.frame(Sample_ID = all_samples, stringsAsFactors = FALSE) %>%
  left_join(meta_info %>% select(Sample_ID, Group), by = "Sample_ID") %>%
  mutate(Group = ifelse(is.na(Group) & grepl("QC", Sample_ID, ignore.case = TRUE), 
                        "QC", Group))
full_meta <- full_meta[match(colnames(expr_final_q), full_meta$Sample_ID), ]

# Perform PCA
pca_input <- t(expr_final_q)
res.pca <- PCA(pca_input, scale.unit = TRUE, graph = FALSE)

# Visualize PCA (PC1 vs PC3)
p_pca <- fviz_pca_ind(
  res.pca, 
  axes = c(1, 3), 
  geom.ind = "point", 
  pointsize = 3,
  col.ind = factor(full_meta$Group, levels = c("Familial", "Sporadic", "QC")), 
  palette = c("#E64B35", "#4DBBD5", "#808080"),
  addEllipses = TRUE,
  ellipse.level = 0.95,
  ellipse.type = "t", 
  legend.title = "Group", 
  mean.point = FALSE
) + 
  theme_minimal() + 
  theme(text = element_text(size = 14), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(title = "PCA after Quantile Normalization")

# Save PCA plot
ggsave(filename = file.path(output_plot, "PCA_normalized_proteomics.png"), 
       plot = p_pca, width = 7, height = 6, dpi = 300)

ggsave(filename = file.path(output_plot, "PCA_normalized_proteomics.svg"), 
       plot = p_pca, width = 7, height = 6, device = "svg")

# ==============================================================================
# 7. Export PCA coordinates
# ==============================================================================
pca_coords <- as.data.frame(res.pca$ind$coord)
pca_coords$Sample_ID <- rownames(pca_coords)
pca_coords$Group <- full_meta$Group

write.csv(pca_coords, 
          file.path(output_dir, "PCA_coordinates.csv"), 
          row.names = FALSE, quote = FALSE)

