# 2026年3月27日 研究日志
今天干了很多事情，以至于过于疲倦了，今天就来简单分享我自己写的用于单细胞测序数据清洗的一些R语言函数吧。


```R
# Check the number of NA and Empty in the meta.table
Check_NA_Empty = function(scdata){
  #require package
  require(readr)
  
  #Get meta data from Seurat object
  meta.table <- scdata@meta.data
  
  #create direction
  if (!dir.exists("./report/NAreport")) dir.create("./report/NAreport", recursive = TRUE)
  sample_name <- deparse(substitute(scdata))
  report_dir <- paste0("./report/NAreport/", sample_name, "_NAreport.csv")
  
  #Count NA and empty value
  na_counts <- colSums(is.na(meta.table))
  empty_counts <- apply(meta.table, 2, function(x){
    x_char = as.character(x)
    sum(x_char == "", na.rm = TRUE)
  })
  
  total_missing <- na_counts + empty_counts
  
  all_na <- na_counts == nrow(meta.table)
  all_empty <- empty_counts == (nrow(meta.table) - na_counts)
  
  #Make final result table for report the NA and empty and save it as CSV file
  result_df <- data.frame(
    Feature = names(na_counts),
    NA_count = na_counts,
    Empty_count = empty_counts,
    Total_missing = total_missing,
    All_NA = all_na,
    All_Empty = all_empty,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  write_csv(result_df, report_dir)
  return(result_df)
}
```

↑**Check_NA_Empty(scdata)** 这是一个用于检查单细胞 Seurat 对象元数据（meta.data）中NA和空值情况的辅助函数。它的主要作用是帮助你快速发现元数据表中是否存在缺失信息，并自动生成报告文件。根据报告你可以看到哪些列完全缺失（通常这意味着该列系统自动生成但实际上并未测量）而哪些列混合着NA或空字符串（这种数据则要根据情况小心判断原因和可能对研究的影响）简单来说这是一个辅助你检查数据集质量的函数。


```R
# Make a Frequency count(One Column Selection)
Frequency_Count = function(scdata, column, sample_name = NULL) {
  #require package
  require(dplyr)
  require(readr)
  
  # Safe column
  safe_column_name <- function(name) {
    gsub("[^[:alnum:]_]", "_", name)
  }
  
  # Get metadata
  meta.table <- scdata@meta.data
  
  # Get column names (raw)
  column_name <- if (is.character(column)) {
    column
  } else {
    deparse(substitute(column))
  }
  
  # Verify that the column exists
  if (!column_name %in% colnames(meta.table)) {
    stop("Column '", column_name, "' not found in metadata")
  }
  
  # Determine the sample name (the passed-in name is used first)
  if (is.null(sample_name)) {
    sample_name <- deparse(substitute(scdata))
  }
  
  # Create a directory (using a safe name)
  safe_sample_name <- safe_column_name(sample_name)
  sample_dir <- file.path("./report", "Frequencyreport", safe_sample_name)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Count using original column name
  counts_table <- meta.table %>%
    dplyr::count(dplyr::across(dplyr::all_of(column_name)), name = "Counts") %>%
    dplyr::mutate(Percentage = round(Counts / sum(Counts) * 100, 2)) %>%
    dplyr::arrange(dplyr::desc(Counts)) %>%
    dplyr::rename(!!column_name := dplyr::all_of(column_name))
  
  # Generate safe file names
  safe_name <- safe_column_name(column_name)
  table_path <- file.path(sample_dir, paste0(safe_name, "_freq.csv"))
  
  # Save CSV
  readr::write_csv(counts_table, table_path)
  
  cat("Frequency table saved to:", table_path, "\n")
  return(counts_table)
}

# Make a Frequency count table list(Manual Selection)
Frequency_Count_List = function(scdata, columns) {
  # Get the name of the object
  sample_name <- deparse(substitute(scdata))
  
  result_list <- lapply(columns, function(col) {
    # Pass the sample name to Frequency_Count
    Frequency_Count(scdata, column = col, sample_name = sample_name)
  })
  
  names(result_list) <- columns
  return(result_list)
}

# Make Frequency count table(Automatic Selection)
Frequency_Count_All = function(scdata, max_unique = 50) {
  # Safe Column
  safe_column_name <- function(name) {
    gsub("[^[:alnum:]_]", "_", name)
  }
  
  # Get the name of the passed object
  sample_name <- deparse(substitute(scdata))
  safe_sample_name <- safe_column_name(sample_name)
  
  # Create Directory Path
  sample_dir <- file.path("./report", "Frequencyreport", safe_sample_name)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  meta.table <- scdata@meta.data
  all_cols <- colnames(meta.table)
  
  result_list <- list()
  skipped_cols <- character(0)
  unique_list <- list()
  
  for (col in all_cols) {
    unique_values <- unique(meta.table[[col]])
    num_unique <- length(unique_values)
    
    if (num_unique > max_unique) {
      skipped_cols <- c(skipped_cols, col)
      next
    }
    
    if (num_unique == 1) {
      unique_list[[col]] <- data.frame(
        Name = safe_column_name(col),
        Value = as.character(unique_values),
        stringsAsFactors = FALSE
      )
    } else {
      tryCatch({
        # Pass the sample name to Frequency_Count
        counts_table <- Frequency_Count(scdata, col, sample_name = sample_name)
        result_list[[safe_column_name(col)]] <- counts_table
      }, error = function(e) {
        warning("Error processing column '", col, "': ", e$message)
        skipped_cols <<- c(skipped_cols, col)
      })
    }
  }
```

**Frequency_Count**家族：包含**Frequency_Count/Frequency_Count_List/Frequency_Count_All**，顾名思义是用来统计频数，生成频数报告的函数。其中**Frequency_Count**是统计单列数值频数的函数，所以需要Seurat数据及列名。**Frequency_Count_List**则是统计你手动选择的列名列表中的所有列的频数。最后**Frequency_Count_All**则是自动化统计所有列（那些ID列等无需统计频数的列用`max_unique`参数进行筛选）要注意的是，**Frequency_Count_List**和**Frequency_Count_All**依赖于底层函数**Frequency_Count**，因此如果需要使用请务必全部粘贴于您的脚本中。

*今天就到这里了，好累，去休息了，下个工作日见咯。*
