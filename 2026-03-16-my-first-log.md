# 2026年3月16日 研究日志
今日主题：解析文献*Complement Factor H (Y402H) Polymorphism Displays Characteristic Features of Age-Related Macular Degeneration and Indicates a Beneficial Role for UV Light Exposure*

今日研究重点：了解该研究的目的及大致实验方法，探究数据类型及数据结构

研究个人摘要：此研究的目的是探究补体因子H多态性Y402H对年龄相关性黄斑变性（AMD）的影响。使用的材料是诱导多能干细胞iPSC，制取方法：从两名无AMD且低风险基因型的受试者和两名患有晚期AMD且高风险基因型的患者中获取皮肤成纤维细胞，通过病毒转导OCT4、SOX2、KLF4和c-MYC这四种转录因子，将其重编程为iPSCs，之后使用RPE分化培养基使其分化为RPE细胞。研究进行了RNA测序分析，并将其发布于GEO数据库，数据集编号GSE91087。

数据可用性：该研究进行了bulk RNA测序，提供了分析结果数据及测序原始数据FASTQ文件，数据集页面： https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91087

分析采用数据备注：\
GSE91087_DE_genes_anno_1.5Fold_June2015.txt -> 作者提供的差异表达分析结果 \
GSE91087_raw_counts_GRCh38.p13_NCBI.tsv -> NCBI根据测序原始数据自动生成的计数矩阵 \
Human.GRCh38.p13.annot.tsv -> NCBI提供的用于对应基因ID和基因symbol的备注表

## 分析用代码存档：

↓导入用于定位数据位置，读取文件以及差异表达分析的R包，除智人数据库org.Hs.eg.db外均为必要


```R
# ----Library----
library(here)
library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)
library(DESeq2)

library(org.Hs.eg.db)
```

↓导入数据并调整数据类型


```R
# ----Data----
authorDEA_txt <- here("Dataset", "GSE91087", "GSE91087_DE_genes_anno_1.5Fold_June2015.txt")
NCBI_Rawcount_tsv <- here("Dataset", "GSE91087", "GSE91087_raw_counts_GRCh38.p13_NCBI.tsv")
annot_tsv <- here("Dataset", "GSE91087", "Human.GRCh38.p13.annot.tsv")

counts <- read_tsv(
  file = NCBI_Rawcount_tsv,
  col_names = TRUE,  
  show_col_types = FALSE 
)

annot <- read_tsv(
  file = annot_tsv,
  col_names = TRUE,
  show_col_types = FALSE
)
annot$GeneID <- as.character(annot$GeneID)
```

↓导入作者差异表达分析结果并进行可视化


```R
#### Author DEA result ####
authorDEA <- read.delim(authorDEA_txt,
                        header = TRUE,
                        sep = "\t",
                        stringsAsFactors = FALSE)

authorDEA$Significance <- ifelse(authorDEA$padj < 0.05 & authorDEA$log2FoldChange > 1, "Up",
                                 ifelse(authorDEA$padj < 0.05 & authorDEA$log2FoldChange < -1, "Down", "Not Sig"))

author_volcano_plot <- ggplot(authorDEA, aes(x = log2FoldChange, y = -log10(padj))) +
  # point(colored by significant)
  geom_point(aes(color = Significance), alpha = 0.6, size = 2) +
  
  # Color set
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "Not Sig" = "gray")) +
  
  # threshold line
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgray") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgray") +
  
  # Gene label(top 10)
  geom_text_repel(
    data = head(authorDEA[order(authorDEA$padj), ], 10),
    aes(label = Gene_Name),
    size = 3, max.overlaps = 15, box.padding = 0.5
  ) +
  
  # title and axis
  labs(
    title = "Differential Expression Analysis",
    x = expression(Log[2] ~ "Fold Change"),
    y = expression(-Log[10] ~ "(Adjusted p-value)"),
    color = "Significance"
  ) +
  
  # theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )
```

↓使用计数矩阵进行差异表达分析（低风险vs高风险）


```R
#### Analyze from count matrix ####
# merge annotate and counts table
counts$GeneID <- rownames(counts)
merged <- counts %>%
  inner_join(annot, by = "GeneID") %>%
  distinct(Symbol, .keep_all = TRUE)

merged <- as.data.frame(merged)
rownames(merged) <- merged$Symbol
merged$Symbol <- NULL
merged$GeneID <- NULL
counts <- merged[,1:6]

# DEA
condition <- factor(c("lowrisk","lowrisk","lowrisk",
                      "highrisk","highrisk","highrisk"))
coldata <- data.frame(row.names = colnames(counts),
                      condition = condition)

coldata$condition <- relevel(coldata$condition, ref = "lowrisk")

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts),
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

head(res)
summary(res)

ordered_res <- res[order(res$padj), ]
res_df <- as.data.frame(ordered_res)

res_df$Significance <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 1, "Up",
                              ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, "Down", "Not Sig"))

Volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  # point(colored by significant)
  geom_point(aes(color = Significance), alpha = 0.6, size = 2) +
  
  # Color set
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "Not Sig" = "gray")) +
  
  # threshold line
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgray") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgray") +
  
  # title and axis
  labs(
    title = "Differential Expression Analysis",
    x = expression(Log[2] ~ "Fold Change"),
    y = expression(-Log[10] ~ "(Adjusted p-value)"),
    color = "Significance"
  ) +
  
  # theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

Volcano_plot
```

# Reference
Hallam, D., Collin, J., Bojic, S., Chichagova, V., Buskin, A., Xu, Y., Lafage, L., Otten, E. G., Anyfantis, G., Mellough, C., Przyborski, S., Alharthi, S., Korolchuk, V., Lotery, A., Saretzki, G., McKibbin, M., Armstrong, L., Steel, D., Kavanagh, D., & Lako, M. (2017). An Induced Pluripotent Stem Cell Patient Specific Model of Complement Factor H (Y402H) Polymorphism Displays Characteristic Features of Age-Related Macular Degeneration and Indicates a Beneficial Role for UV Light Exposure. Stem Cells, 35(11), 2305–2320. https://doi.org/10.1002/STEM.2708

*备注：这篇博客是我第一次尝试生成技术博客，说来惭愧，之前从来没写过博客，因此也不知道应当注意些什么。总之计划是每个工作日都要更新研究日志。*
