# 2026年3月17日 研究日志
今日主题：解析文献*Identification of Ferroptosis-Related Gene in Age-Related Macular Degeneration Using Machine Learning*

今日研究重点：了解该研究的目的及大致实验方法，探究数据类型及数据结构

研究个人笔记：此研究使用了GEO数据集GSE29801，这是一个2012年生成的微阵列数据集。说实话数据有些老旧（谁还用微阵列啊），不过胜在该文献对操作步骤描述详细，对如何分析差异表达基因和如何使用机器学习算法都进行了比较详细的描述（当然想要直接复现还是需要更多细节）。对这篇文献进行复现的主要目的是了解差异表达分析到机器学习生成模型的工作流（今日只进行到差异表达分析）。

数据可用性：该研究是使用GSE29801的下游研究。

GSE29801 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29801

## 分析用代码存档：

↓所有可能用得上的库


```R
# Library
library(GEOquery)
library(data.table)
library(limma)
library(dplyr)
library(clusterProfiler)
library(glmnet)
library(e1071)
library(randomForest)
library(pROC)
library(here)
library(ggplot2)
library(ggrepel)
```

↓所有数据，文章实际上使用了数据集中所有黄斑区的样本（毕竟是AMD相关），筛选标准在文献中，这里将其粘贴存入了txt文件中。除了主要数据集，实际上文章还利用了其他三个数据集进行了验证，今天的分析未涉及这一部分。


```R
# ----Data----
# Main dataset: GSE29801
gse_data <- getGEO(GEO = "GSE29801")
gse_data <- gse_data[[1]]

# Validation dataset
validation_set_1 <- getGEO(GEO = "GSE99248")
validation_set_2 <- getGEO(GEO = "GSE50195")
validation_set_3 <- getGEO(GEO = "GSE125564") 

chosen_sample_table <- fread(here("Dataset","GSE29801", "table1.txt"), header = TRUE, stringsAsFactors = FALSE)
chosen_samples <- as.character(chosen_sample_table$`Sample GEO accession`)

Genetic_risk_table <- fread(here("Dataset","GSE29801", "table2.txt"), header = TRUE, stringsAsFactors = FALSE)
```

接下来的分为两个方法分析数据，因为原文献中没有提及他们究竟是先筛选数据再做数据清理，还是先清理再做筛选。因为limma本质是线性模型，样本的不同一定会影响结果，实际上文章也没详细描述预处理参数，因此两个方法都不能得到原版的452个差异基因的结果。（已致信文章通信作者，暂未取得回复）

↓方案1：先筛选再清理数据


```R
# ----Replicate experiment----
# Subset before clean data & limma
all_samples <- rownames(pData(gse_data)) 
selected_samples <- all_samples %in% chosen_samples 
gse_subset <- gse_data[, selected_samples]

print(paste("Original amount:", ncol(gse_data)))
print(paste("After chosen:", ncol(gse_subset)))

### DEGA Step1：Prepare express matrix ####
# Clean micro array data
feature_data <- fData(gse_subset)
valid_genes <- feature_data[!is.na(feature_data$GENE), c("ID", "GENE_SYMBOL")]
expr_matrix <- exprs(gse_subset)

gene_mapping <- setNames(valid_genes$GENE_SYMBOL, valid_genes$ID)
all_probes <- rownames(expr_matrix)
valid_probes <- all_probes[all_probes %in% names(gene_mapping)]

expr_preprocessed <- backgroundCorrect(expr_matrix, method = "normexp", offset = 16) %>%
  log2() %>%
  normalizeBetweenArrays(method = "quantile")

expr_preprocessed_filtered <- expr_preprocessed[valid_probes, ]
probe_symbols <- gene_mapping[rownames(expr_preprocessed_filtered)]

expr_gene_final <- do.call(rbind, lapply(unique(probe_symbols), function(gene) {
  probe_indices <- which(probe_symbols == gene)
  apply(expr_preprocessed_filtered[probe_indices, , drop = FALSE], 2, 
        function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE))
}))

rownames(expr_gene_final) <- unique(probe_symbols)

### DEGA Step2: Process express matrix with limma####
# Group Information
sample_info <- pData(gse_subset)
groups <- sample_info$source_name_ch1 

# design matrix
design <- model.matrix(~0 + groups)
colnames(design) <- paste0("Tissue_", gsub("groups|[- ,]", "", colnames(design)))

# Linear model fitting
fit <- lmFit(expr_gene_final, design) 

# Contrast matrix
contrast_matrix <- makeContrasts(
  AMD_vs_Normal = Tissue_MacularRPEchoroidAMD - Tissue_MacularRPEchoroidnormal,
  levels = design
)

fit_contrasts <- contrasts.fit(fit, contrast_matrix)

fit_eb <- eBayes(fit_contrasts)

# result
results <- topTable(fit_eb, 
                    coef = "AMD_vs_Normal", 
                    number = Inf,  
                    adjust.method = "BH")  

results$gene <- rownames(results)

### DEGA Step3：Analyse DEA result and Draw Volcano plot ####
# Prepare plot data
p_threshold <- 0.05
logFC_threshold <- 0.263

volcano_data <- results %>%
  mutate(
    log10_pval = -log10(P.Value),
    # significant mark
    Significance = case_when(
      P.Value < p_threshold & logFC > logFC_threshold ~ "Up-regulated",
      P.Value < p_threshold & logFC < -logFC_threshold ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

# Check the number of significant genes
up_genes <- sum(volcano_data$Significance == "Up-regulated", na.rm = TRUE)
down_genes <- sum(volcano_data$Significance == "Down-regulated", na.rm = TRUE)

# Up regulate genes top5
top_up <- volcano_data %>%
  filter(Significance == "Up-regulated" & !is.na(gene)) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 5)

# Down regulate genes top5
top_down <- volcano_data %>%
  filter(Significance == "Down-regulated" & !is.na(gene)) %>%
  arrange(logFC) %>%
  slice_head(n = 5)

genes_to_label <- bind_rows(top_up, top_down)

# Main Plot
volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = log10_pval)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
  # Color
  scale_color_manual(
    values = c(
      "Up-regulated" = "red",
      "Down-regulated" = "blue",
      "Not significant" = "grey"
    ),
    breaks = c("Up-regulated", "Down-regulated")
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  
  labs(
    x = expression("log"[2] * " Fold Change"),
    y = expression("-log"[10] * " (P-value)"),
    title = paste0("Volcano Plot of Differential Expression\n",
                   up_genes, " up-regulated, ", down_genes, " down-regulated genes\n",
                   "(p < 0.05, |log2FC| > 0.263)")
  ) +
  
  # threshold line
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(p_threshold), 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  
  # Annotation
  geom_text_repel(
    data = genes_to_label,
    aes(label = gene),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.5,
    segment.color = "grey40",
    segment.size = 0.3,
    min.segment.length = 0.1,
    max.overlaps = 20
  ) +
  
  # threshold label
  annotate("text", x = -logFC_threshold - 0.3, y = max(volcano_data$log10_pval)*0.95, 
           label = "log2FC = -0.263", color = "black", size = 3.3) +
  annotate("text", x = logFC_threshold + 0.3, y = max(volcano_data$log10_pval)*0.95, 
           label = "log2FC = 0.263", color = "black", size = 3.3) +
  annotate("text", x = max(volcano_data$logFC)*0.9, y = -log10(p_threshold) + 0.1, 
           label = "p = 0.05", color = "black", size = 3.3)

volcano_plot
```

↓方案2：先清洗数据后进行筛选


```R
# ----Analysis----
### DEGA Step1：Prepare express matrix ####
# Clean micro array data
feature_data <- fData(gse_data)
valid_genes <- feature_data[!is.na(feature_data$GENE), c("ID", "GENE_SYMBOL")]
expr_matrix <- exprs(gse_data)

gene_mapping <- setNames(valid_genes$GENE_SYMBOL, valid_genes$ID)
all_probes <- rownames(expr_matrix)
valid_probes <- all_probes[all_probes %in% names(gene_mapping)]

expr_preprocessed <- backgroundCorrect(expr_matrix, method = "normexp", offset = 16) %>%
  log2() %>%
  normalizeBetweenArrays(method = "quantile")

# Subset after clean data
expr_preprocessed_filtered <- expr_preprocessed[valid_probes, ]
probe_symbols <- gene_mapping[rownames(expr_preprocessed_filtered)]

expr_gene_final <- do.call(rbind, lapply(unique(probe_symbols), function(gene) {
  probe_indices <- which(probe_symbols == gene)
  apply(expr_preprocessed_filtered[probe_indices, , drop = FALSE], 2, 
        function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE))
}))

rownames(expr_gene_final) <- unique(probe_symbols)
expr_gene_final <- expr_gene_final[, chosen_samples, drop = FALSE]

### DEGA Step2: Process express matrix with limma####
# Group Information
sample_info <- pData(gse_subset)
groups <- sample_info$source_name_ch1 

# design matrix
design <- model.matrix(~0 + groups)
colnames(design) <- paste0("Tissue_", gsub("groups|[- ,]", "", colnames(design)))

# Linear model fitting
fit <- lmFit(expr_gene_final, design) 

# Contrast matrix
contrast_matrix <- makeContrasts(
  AMD_vs_Normal = Tissue_MacularRPEchoroidAMD - Tissue_MacularRPEchoroidnormal,
  levels = design
)

fit_contrasts <- contrasts.fit(fit, contrast_matrix)

fit_eb <- eBayes(fit_contrasts)

# result
results <- topTable(fit_eb, 
                    coef = "AMD_vs_Normal", 
                    number = Inf,  
                    adjust.method = "BH")  

results$gene <- rownames(results)

### DEGA Step3：Analyse DEA result and Draw Volcano plot ####
# Prepare plot data
p_threshold <- 0.05
logFC_threshold <- 0.263

volcano_data <- results %>%
  mutate(
    log10_pval = -log10(P.Value),
    # significant mark
    Significance = case_when(
      P.Value < p_threshold & logFC > logFC_threshold ~ "Up-regulated",
      P.Value < p_threshold & logFC < -logFC_threshold ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

# Check the number of significant genes 
up_genes <- sum(volcano_data$Significance == "Up-regulated", na.rm = TRUE)
down_genes <- sum(volcano_data$Significance == "Down-regulated", na.rm = TRUE)

# Up regulate genes top5
top_up <- volcano_data %>%
  filter(Significance == "Up-regulated" & !is.na(gene)) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 5)

# Down regulate genes top5
top_down <- volcano_data %>%
  filter(Significance == "Down-regulated" & !is.na(gene)) %>%
  arrange(logFC) %>%
  slice_head(n = 5)

genes_to_label <- bind_rows(top_up, top_down)

# Main Plot
volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = log10_pval)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
  # Color
  scale_color_manual(
    values = c(
      "Up-regulated" = "red",
      "Down-regulated" = "blue",
      "Not significant" = "grey"
    ),
    breaks = c("Up-regulated", "Down-regulated")
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  
  labs(
    x = expression("log"[2] * " Fold Change"),
    y = expression("-log"[10] * " (P-value)"),
    title = paste0("Volcano Plot of Differential Expression\n",
                   up_genes, " up-regulated, ", down_genes, " down-regulated genes\n",
                   "(p < 0.05, |log2FC| > 0.263)")
  ) +
  
  # threshold line
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(p_threshold), 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  
  # Annotation
  geom_text_repel(
    data = genes_to_label,
    aes(label = gene),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.5,
    segment.color = "grey40",
    segment.size = 0.3,
    min.segment.length = 0.1,
    max.overlaps = 20
  ) +
  
  # threshold label 
  annotate("text", x = -logFC_threshold - 0.3, y = max(volcano_data$log10_pval)*0.95, 
           label = "log2FC = -0.263", color = "black", size = 3.3) +
  annotate("text", x = logFC_threshold + 0.3, y = max(volcano_data$log10_pval)*0.95, 
           label = "log2FC = 0.263", color = "black", size = 3.3) +
  annotate("text", x = max(volcano_data$logFC)*0.9, y = -log10(p_threshold) + 0.1, 
           label = "p = 0.05", color = "black", size = 3.3)

volcano_plot
```

# Reference
Zhu, M., & Yu, J. (2024). Identification of Ferroptosis-Related Gene in Age-Related Macular Degeneration Using Machine Learning. Immunity, Inflammation and Disease , 12(12). https://doi.org/10.1002/iid3.70059

*备注：处理微阵列数据真是费心费神，今天先到这里了，后续的机器学习部分等作者回信（不知道等的等不到）再说吧。明天想先总结一下堆积如山的未整理的生物信息学相关的python函数*
