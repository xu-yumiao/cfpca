# ==============================================================================
# GO Enrichment Analysis for WGCNA Modules (BP/CC/MF)
# Description: Perform three-dimensional GO enrichment analysis (Biological 
#              Process, Cellular Component, Molecular Function) for selected
#              WGCNA modules
# ==============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# ==============================================================================
# 1. Define paths and load data
# ==============================================================================
output_data <- "results/proteomics/WGCNA"

load(file.path(output_data, "step5_module_analysis.RData"))

# ==============================================================================
# 2. Define modules for enrichment analysis
# ==============================================================================
# Specify modules of interest (replace with selected module colors)
key_modules <- c("module1", "module2", "module3")

# ==============================================================================
# 3. Set enrichment analysis parameters
# ==============================================================================
pvalue_cutoff <- 0.05
qvalue_cutoff <- 0.2

# Initialize result storage
GO_BP_results <- list()
GO_CC_results <- list()
GO_MF_results <- list()

# ==============================================================================
# 4. Perform three-dimensional GO enrichment for each module
# ==============================================================================
for(module_color in key_modules) {
  
  # Extract module proteins
  module_proteins <- module_analysis_results[[module_color]]$Protein
  
  # GO-BP enrichment
  tryCatch({
    ego_bp <- enrichGO(gene = module_proteins,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = pvalue_cutoff,
                       qvalueCutoff = qvalue_cutoff,
                       readable = TRUE)
    
    if(!is.null(ego_bp) && nrow(ego_bp@result) > 0) {
      GO_BP_results[[module_color]] <- ego_bp
      
      go_bp_df <- as.data.frame(ego_bp)
      write.csv(go_bp_df,
                file.path(output_data, sprintf("module_%s_GO_BP_full.csv", module_color)),
                row.names = FALSE)
    } else {
      GO_BP_results[[module_color]] <- NULL
    }
  }, error = function(e) {
    GO_BP_results[[module_color]] <- NULL
  })
  
  # GO-CC enrichment
  tryCatch({
    ego_cc <- enrichGO(gene = module_proteins,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff = pvalue_cutoff,
                       qvalueCutoff = qvalue_cutoff,
                       readable = TRUE)
    
    if(!is.null(ego_cc) && nrow(ego_cc@result) > 0) {
      GO_CC_results[[module_color]] <- ego_cc
      
      go_cc_df <- as.data.frame(ego_cc)
      write.csv(go_cc_df,
                file.path(output_data, sprintf("module_%s_GO_CC_full.csv", module_color)),
                row.names = FALSE)
    } else {
      GO_CC_results[[module_color]] <- NULL
    }
  }, error = function(e) {
    GO_CC_results[[module_color]] <- NULL
  })
  
  # GO-MF enrichment
  tryCatch({
    ego_mf <- enrichGO(gene = module_proteins,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff = pvalue_cutoff,
                       qvalueCutoff = qvalue_cutoff,
                       readable = TRUE)
    
    if(!is.null(ego_mf) && nrow(ego_mf@result) > 0) {
      GO_MF_results[[module_color]] <- ego_mf
      
      go_mf_df <- as.data.frame(ego_mf)
      write.csv(go_mf_df,
                file.path(output_data, sprintf("module_%s_GO_MF_full.csv", module_color)),
                row.names = FALSE)
    } else {
      GO_MF_results[[module_color]] <- NULL
    }
  }, error = function(e) {
    GO_MF_results[[module_color]] <- NULL
  })
}

# ==============================================================================
# 5. Save enrichment results
# ==============================================================================
save(GO_BP_results, GO_CC_results, GO_MF_results,
     key_modules,
     file = file.path(output_data, "GO_three_dimensions_results.RData"))

# ==============================================================================
# 6. Generate enrichment summary statistics
# ==============================================================================
summary_stats <- data.frame()

for(module_color in key_modules) {
  
  n_bp <- ifelse(is.null(GO_BP_results[[module_color]]), 0,
                 nrow(GO_BP_results[[module_color]]@result))
  n_cc <- ifelse(is.null(GO_CC_results[[module_color]]), 0,
                 nrow(GO_CC_results[[module_color]]@result))
  n_mf <- ifelse(is.null(GO_MF_results[[module_color]]), 0,
                 nrow(GO_MF_results[[module_color]]@result))
  
  summary_stats <- rbind(summary_stats, data.frame(
    Module = module_color,
    Protein_Count = nrow(module_analysis_results[[module_color]]),
    GO_BP_Terms = n_bp,
    GO_CC_Terms = n_cc,
    GO_MF_Terms = n_mf,
    Total_GO_Terms = n_bp + n_cc + n_mf
  ))
}

write.csv(summary_stats,
          file.path(output_data, "GO_enrichment_summary.csv"),
          row.names = FALSE)

# ==============================================================================
# End of script
# ==============================================================================
