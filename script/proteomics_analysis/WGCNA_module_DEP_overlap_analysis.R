# ============================================================
# WGCNA 第六步：模块富集分析与差异表达交叉验证
# ============================================================

library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(svglite)
library(enrichplot)

# 加载数据
output_data <- "C:/Users/chevy/Desktop/科研/前列腺癌文章修改/处理数据/蛋白组/WGCNA"
output_plot <- "C:/Users/chevy/Desktop/科研/前列腺癌文章修改/原始图片/WGCNA"
input_limma <- "C:/Users/chevy/Desktop/科研/前列腺癌文章修改/处理数据/蛋白组/Limma_Diff_Results_Ultimate.csv"

load(file.path(output_data, "step5_module_analysis.RData"))

cat("数据已加载\n")

# ============================================================
# 1. 读取差异表达结果
# ============================================================

limma_results <- read.csv(input_limma)
cat("\n差异表达结果:\n")
cat("总蛋白数:", nrow(limma_results), "\n")
cat("上调蛋白 (logFC>0, adj.P<0.05):", sum(limma_results$adj.P.Val < 0.05 & limma_results$logFC > 0), "\n")
cat("下调蛋白 (logFC<0, adj.P<0.05):", sum(limma_results$adj.P.Val < 0.05 & limma_results$logFC < 0), "\n")

# 提取显著差异蛋白
sig_threshold <- 0.05
DEPs <- limma_results %>%
  filter(adj.P.Val < sig_threshold) %>%
  mutate(Direction = ifelse(logFC > 0, "Up", "Down"))

cat("显著差异蛋白总数:", nrow(DEPs), "\n")

# ============================================================
# 2. 模块与差异表达的交叉验证
# ============================================================

cat("\n========== 模块与差异表达交叉统计 ==========\n")

significant_modules_list <- c("turquoise", "black", "greenyellow", "blue", "brown")

overlap_summary <- data.frame()

for(module_color in significant_modules_list) {
  
  # 提取模块蛋白
  module_proteins <- module_analysis_results[[module_color]]$Protein
  n_module <- length(module_proteins)
  
  # 与差异蛋白交叉
  overlap_DEPs <- DEPs[DEPs$Protein %in% module_proteins, ]
  n_overlap <- nrow(overlap_DEPs)
  n_up <- sum(overlap_DEPs$Direction == "Up")
  n_down <- sum(overlap_DEPs$Direction == "Down")
  
  # 超几何检验：模块是否富集差异蛋白
  total_proteins <- length(unique(mergedColors))  # 所有蛋白数
  total_DEPs <- nrow(DEPs)  # 总差异蛋白数
  
  # 使用phyper进行超几何检验
  # phyper(q, m, n, k, lower.tail=FALSE)
  # q: 模块中差异蛋白数-1
  # m: 总差异蛋白数
  # n: 总非差异蛋白数
  # k: 模块蛋白数
  enrich_pval <- phyper(n_overlap - 1, 
                        total_DEPs, 
                        total_proteins - total_DEPs, 
                        n_module, 
                        lower.tail = FALSE)
  
  overlap_summary <- rbind(overlap_summary, data.frame(
    Module = module_color,
    Module_Size = n_module,
    Overlap_DEPs = n_overlap,
    Overlap_Up = n_up,
    Overlap_Down = n_down,
    Overlap_Percent = sprintf("%.1f%%", 100*n_overlap/n_module),
    Enrichment_Pvalue = enrich_pval
  ))
  
  cat(sprintf("\n%s模块:\n", module_color))
  cat(sprintf("  模块蛋白数: %d\n", n_module))
  cat(sprintf("  重叠差异蛋白: %d (%.1f%%)\n", n_overlap, 100*n_overlap/n_module))
  cat(sprintf("    上调: %d\n", n_up))
  cat(sprintf("    下调: %d\n", n_down))
  cat(sprintf("  富集P值: %.2e\n", enrich_pval))
}

# 保存交叉统计结果
write.csv(overlap_summary,
          file.path(output_data, "module_DEP_overlap_summary.csv"),
          row.names = FALSE)

cat("\n交叉统计结果已保存\n")

# ============================================================
# 3. 准备富集分析
# ============================================================

cat("\n========== 开始富集分析 ==========\n")

# 设置富集分析参数
pvalue_cutoff <- 0.05
qvalue_cutoff <- 0.2

# 初始化结果列表
GO_results <- list()
KEGG_results <- list()

# ============================================================
# 4. 对每个模块进行GO-BP富集
# ============================================================

for(module_color in significant_modules_list) {
  
  cat(sprintf("\n处理模块: %s\n", module_color))
  
  # 提取模块蛋白名称
  module_proteins <- module_analysis_results[[module_color]]$Protein
  
  cat(sprintf("  模块蛋白数: %d\n", length(module_proteins)))
  
  # GO-BP富集分析
  cat("  进行GO-BP富集分析...\n")
  
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
      GO_results[[module_color]] <- ego_bp
      cat(sprintf("    富集到%d个GO-BP条目\n", nrow(ego_bp@result)))
    } else {
      cat("    未富集到显著GO-BP条目\n")
      GO_results[[module_color]] <- NULL
    }
  }, error = function(e) {
    cat(sprintf("    GO-BP富集失败: %s\n", e$message))
    GO_results[[module_color]] <- NULL
  })
  
  # KEGG通路富集分析
  cat("  进行KEGG富集分析...\n")
  
  tryCatch({
    # 将基因符号转换为Entrez ID
    gene_entrez <- bitr(module_proteins, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db)
    
    cat(sprintf("    成功转换%d个蛋白为Entrez ID\n", nrow(gene_entrez)))
    
    if(nrow(gene_entrez) > 10) {  # 至少需要10个基因
      ekegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                          organism = "hsa",
                          pvalueCutoff = pvalue_cutoff,
                          qvalueCutoff = qvalue_cutoff)
      
      if(!is.null(ekegg) && nrow(ekegg@result) > 0) {
        KEGG_results[[module_color]] <- ekegg
        cat(sprintf("    富集到%d个KEGG通路\n", nrow(ekegg@result)))
      } else {
        cat("    未富集到显著KEGG通路\n")
        KEGG_results[[module_color]] <- NULL
      }
    } else {
      cat("    转换的基因数太少，跳过KEGG分析\n")
      KEGG_results[[module_color]] <- NULL
    }
  }, error = function(e) {
    cat(sprintf("    KEGG富集失败: %s\n", e$message))
    KEGG_results[[module_color]] <- NULL
  })
}

# ============================================================
# 5. 保存富集结果表格
# ============================================================

cat("\n保存富集结果...\n")

for(module_color in significant_modules_list) {
  
  # 保存GO结果
  if(!is.null(GO_results[[module_color]])) {
    go_df <- as.data.frame(GO_results[[module_color]])
    write.csv(go_df,
              file.path(output_data, sprintf("module_%s_GO_BP.csv", module_color)),
              row.names = FALSE)
  }
  
  # 保存KEGG结果
  if(!is.null(KEGG_results[[module_color]])) {
    kegg_df <- as.data.frame(KEGG_results[[module_color]])
    write.csv(kegg_df,
              file.path(output_data, sprintf("module_%s_KEGG.csv", module_color)),
              row.names = FALSE)
  }
}

cat("富集结果表格已保存\n")

# ============================================================
# 6. 可视化：GO-BP气泡图
# ============================================================

cat("\n绘制GO-BP富集气泡图...\n")

# 绘图参数设置区
plot_width <- 14
plot_height <- 10
base_size <- 11
title_size <- 13
axis_title_size <- 11
axis_text_size <- 9
point_size_range <- c(3, 8)

# 合并所有模块的GO结果
go_plot_data <- data.frame()

for(module_color in significant_modules_list) {
  if(!is.null(GO_results[[module_color]])) {
    go_df <- as.data.frame(GO_results[[module_color]])
    # 取每个模块top 10
    go_df_top <- go_df %>%
      arrange(p.adjust) %>%
      head(10) %>%
      mutate(Module = module_color)
    go_plot_data <- rbind(go_plot_data, go_df_top)
  }
}

if(nrow(go_plot_data) > 0) {
  
  # 截短过长的描述
  go_plot_data$Description_short <- ifelse(nchar(go_plot_data$Description) > 60,
                                           paste0(substr(go_plot_data$Description, 1, 57), "..."),
                                           go_plot_data$Description)
  
  # 按模块和P值排序
  go_plot_data <- go_plot_data %>%
    arrange(Module, p.adjust) %>%
    mutate(Description_short = factor(Description_short, levels = unique(Description_short)))
  
  p_go <- ggplot(go_plot_data, aes(x = Module, y = Description_short)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "#B2182B", high = "#2166AC",
                         name = "Adjusted\nP value") +
    scale_size_continuous(range = point_size_range, name = "Gene\nCount") +
    labs(x = "", y = "", title = "GO Biological Process Enrichment") +
    theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = axis_text_size, angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = axis_text_size),
      axis.title = element_text(size = axis_title_size),
      legend.title = element_text(size = axis_text_size, face = "bold"),
      legend.text = element_text(size = axis_text_size - 1),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_plot, "06_GO_BP_bubble_plot.png"),
         p_go, width = plot_width, height = plot_height, dpi = 600)
  
  ggsave(file.path(output_plot, "06_GO_BP_bubble_plot.svg"),
         p_go, width = plot_width, height = plot_height)
  
  cat("GO-BP气泡图已保存\n")
} else {
  cat("警告，没有GO富集结果可绘图\n")
}

# ============================================================
# 7. 可视化：KEGG气泡图
# ============================================================

cat("\n绘制KEGG富集气泡图...\n")

kegg_plot_data <- data.frame()

for(module_color in significant_modules_list) {
  if(!is.null(KEGG_results[[module_color]])) {
    kegg_df <- as.data.frame(KEGG_results[[module_color]])
    kegg_df_top <- kegg_df %>%
      arrange(p.adjust) %>%
      head(10) %>%
      mutate(Module = module_color)
    kegg_plot_data <- rbind(kegg_plot_data, kegg_df_top)
  }
}

if(nrow(kegg_plot_data) > 0) {
  
  # 截短过长的描述
  kegg_plot_data$Description_short <- ifelse(nchar(kegg_plot_data$Description) > 60,
                                             paste0(substr(kegg_plot_data$Description, 1, 57), "..."),
                                             kegg_plot_data$Description)
  
  kegg_plot_data <- kegg_plot_data %>%
    arrange(Module, p.adjust) %>%
    mutate(Description_short = factor(Description_short, levels = unique(Description_short)))
  
  p_kegg <- ggplot(kegg_plot_data, aes(x = Module, y = Description_short)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "#B2182B", high = "#2166AC",
                         name = "Adjusted\nP value") +
    scale_size_continuous(range = point_size_range, name = "Gene\nCount") +
    labs(x = "", y = "", title = "KEGG Pathway Enrichment") +
    theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = axis_text_size, angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = axis_text_size),
      axis.title = element_text(size = axis_title_size),
      legend.title = element_text(size = axis_text_size, face = "bold"),
      legend.text = element_text(size = axis_text_size - 1),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_plot, "06_KEGG_bubble_plot.png"),
         p_kegg, width = plot_width, height = plot_height, dpi = 600)
  
  ggsave(file.path(output_plot, "06_KEGG_bubble_plot.svg"),
         p_kegg, width = plot_width, height = plot_height)
  
  cat("KEGG气泡图已保存\n")
} else {
  cat("警告：没有KEGG富集结果可绘图\n")
}

# ============================================================
# 8. 保存R数据对象
# ============================================================

save(GO_results, KEGG_results, overlap_summary, DEPs,
     file = file.path(output_data, "step6_enrichment_analysis.RData"))

cat("\n中间数据已保存\n")

# ============================================================
# 9. 输出摘要
# ============================================================

cat("\n========== 第六步完成摘要 ==========\n")
cat("富集分析参数:\n")
cat("  P值阈值:", pvalue_cutoff, "\n")
cat("  FDR阈值:", qvalue_cutoff, "\n\n")

cat("各模块富集结果:\n")
for(module_color in significant_modules_list) {
  n_go <- ifelse(is.null(GO_results[[module_color]]), 0, nrow(GO_results[[module_color]]@result))
  n_kegg <- ifelse(is.null(KEGG_results[[module_color]]), 0, nrow(KEGG_results[[module_color]]@result))
  cat(sprintf("  %s: GO-BP=%d, KEGG=%d\n", module_color, n_go, n_kegg))
}

cat("\n生成的文件:\n")
cat("  - module_DEP_overlap_summary.csv: 模块与差异表达交叉统计\n")
cat("  - module_<color>_GO_BP.csv: 各模块GO-BP富集结果\n")
cat("  - module_<color>_KEGG.csv: 各模块KEGG富集结果\n")
cat("  - 06_GO_BP_bubble_plot.png/svg: GO-BP气泡图\n")
cat("  - 06_KEGG_bubble_plot.png/svg: KEGG气泡图\n")
  
  
  
  