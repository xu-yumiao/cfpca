# ==============================================================================
# GO-BP Enrichment Analysis for Differentially Expressed Proteins
# Description: Perform GO Biological Process enrichment analysis separately 
#              for upregulated and downregulated DEPs and visualize top 15 
#              enriched pathways
# ==============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(svglite)

set.seed(123)

# ==============================================================================
# 1. Define paths
# ==============================================================================
DEP_file <- "results/proteomics/differential_proteins_results.csv"
output_dir <- "results/proteomics/enrichment"
output_plot <- "figures/enrichment"

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_plot, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Load DEPs and separate by direction
# ==============================================================================
DEP_data <- read.csv(DEP_file, stringsAsFactors = FALSE)

# Separate upregulated and downregulated DEPs
DEP_up <- DEP_data %>%
  filter(adj.P.Val < 0.05 & logFC > log2(1.5)) %>%
  pull(Protein)

DEP_down <- DEP_data %>%
  filter(adj.P.Val < 0.05 & logFC < -log2(1.5)) %>%
  pull(Protein)

# ==============================================================================
# 3. Convert gene symbols to Entrez IDs
# ==============================================================================
convert_to_entrez <- function(gene_symbols) {
  entrez_ids <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  return(entrez_ids)
}

DEP_up_entrez <- convert_to_entrez(DEP_up)
DEP_down_entrez <- convert_to_entrez(DEP_down)

# ==============================================================================
# 4. GO-BP enrichment analysis
# ==============================================================================
# Upregulated DEPs
GO_BP_up <- enrichGO(
  gene = DEP_up_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Downregulated DEPs
GO_BP_down <- enrichGO(
  gene = DEP_down_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ==============================================================================
# 5. Save enrichment results
# ==============================================================================
write.csv(as.data.frame(GO_BP_up), 
          file.path(output_dir, "GO_BP_upregulated_DEPs.csv"), 
          row.names = FALSE)

write.csv(as.data.frame(GO_BP_down), 
          file.path(output_dir, "GO_BP_downregulated_DEPs.csv"), 
          row.names = FALSE)

# ==============================================================================
# 6. Visualization function for enrichment bubble plot
# ==============================================================================
plot_enrichment_bubble <- function(enrich_result, title, top_n = 15) {
  
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    return(NULL)
  }
  
  plot_data <- as.data.frame(enrich_result) %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    mutate(
      GeneRatio_numeric = sapply(GeneRatio, function(x) eval(parse(text = x))),
      Description = factor(Description, levels = rev(Description))
    )
  
  p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(
      low = "#A50F15",
      high = "#FEE5D9",
      name = "Adjusted\nP-value"
    ) +
    scale_size_continuous(
      range = c(3, 8),
      name = "Gene\nCount"
    ) +
    labs(
      title = title,
      x = "Gene Ratio",
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
    )
  
  return(p)
}

# ==============================================================================
# 7. Generate enrichment plots
# ==============================================================================
p_GO_BP_up <- plot_enrichment_bubble(
  GO_BP_up,
  title = "GO-BP Enrichment: Upregulated DEPs"
)

p_GO_BP_down <- plot_enrichment_bubble(
  GO_BP_down,
  title = "GO-BP Enrichment: Downregulated DEPs"
)

# ==============================================================================
# 8. Save plots
# ==============================================================================
save_plot <- function(plot_obj, filename_base, width = 10, height = 8, dpi = 600) {
  if (is.null(plot_obj)) {
    return()
  }
  
  # PNG
  ggsave(
    filename = file.path(output_plot, paste0(filename_base, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  # SVG
  ggsave(
    filename = file.path(output_plot, paste0(filename_base, ".svg")),
    plot = plot_obj,
    width = width,
    height = height,
    device = svglite,
    bg = "white"
  )
}

save_plot(p_GO_BP_up, "GO_BP_upregulated_DEPs_top15")
save_plot(p_GO_BP_down, "GO_BP_downregulated_DEPs_top15")

# ==============================================================================
# End of script
# ==============================================================================
 