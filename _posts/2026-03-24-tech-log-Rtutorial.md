# 2026年3月24日 研究日志
数据清洗中，具体内容需保密，记录一下清洗用的函数。说真的，学习的时候都用的预处理的干净数据集，真到了现实环境中才发现，什么稀奇古怪的问题都有。干净的数据才少见，脏的才到处都是，所以逢山开路遇水架桥的本事还是得有，这里介绍一点R语言中的工具。

↓这些库得下，我猜你应该有下吧。总之先放在这里了。


```R
library(janitor)
library(tidyr)
library(dplyr)
library(stringr)
```

    
    Attaching package: 'janitor'
    
    
    The following objects are masked from 'package:stats':
    
        chisq.test, fisher.test
    
    
    
    Attaching package: 'dplyr'
    
    
    The following objects are masked from 'package:stats':
    
        filter, lag
    
    
    The following objects are masked from 'package:base':
    
        intersect, setdiff, setequal, union
    
    
    

## 统一列名的方法: *janitor::clean_names()*

我们现在有一个原始数据集df_raw, 它长这个样子：


```R
df_raw <- data.frame(
  "Sample ID" = 1:3,
  "Gene.Name" = c("TP53", "BRCA1", "EGFR"),
  "Expression Level (%)" = c(12.5, 8.3, 15.1),
  "  Tissue Type " = c("Lung", "Breast", "Colon")
)
```


```R
df_raw
```


<table class="dataframe">
<caption>A data.frame: 3 × 4</caption>
<thead>
	<tr><th scope=col>Sample.ID</th><th scope=col>Gene.Name</th><th scope=col>Expression.Level....</th><th scope=col>X..Tissue.Type.</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>1</td><td>TP53 </td><td>12.5</td><td>Lung  </td></tr>
	<tr><td>2</td><td>BRCA1</td><td> 8.3</td><td>Breast</td></tr>
	<tr><td>3</td><td>EGFR </td><td>15.1</td><td>Colon </td></tr>
</tbody>
</table>



可以看到它们有的用空格作为分隔，有的用点作为分隔，有的甚至在开头有莫名其妙的空格。这并非不会发生，在多人参与数据收集的大型数据集中，这种现象经常发生。可怖的空格可能让你陷入到不明原因的bug中。由于问题出现在数据集中，你怎么检查自己的代码都不会找到解决办法的。因此，在开始对数据集进行分析前，我们可以进行以下操作：


```R
df_clean <- df_raw |> janitor::clean_names()
```


```R
df_clean
```


<table class="dataframe">
<caption>A data.frame: 3 × 4</caption>
<thead>
	<tr><th scope=col>sample_id</th><th scope=col>gene_name</th><th scope=col>expression_level</th><th scope=col>x_tissue_type</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>1</td><td>TP53 </td><td>12.5</td><td>Lung  </td></tr>
	<tr><td>2</td><td>BRCA1</td><td> 8.3</td><td>Breast</td></tr>
	<tr><td>3</td><td>EGFR </td><td>15.1</td><td>Colon </td></tr>
</tbody>
</table>



snack_case \o/\o/\o/\o/

## 找出未匹配样本的方法: *dplyr::anti_join()*

来看这两个数据集，它们一个代表了元数据集，另一个代表了表达量数据集


```R
df_meta <- data.frame(
  sample_id = c("S1", "S2", "S3", "S4"),
  tissue = c("Lung", "Breast", "Colon", "Liver")
)

df_expr <- data.frame(
  sample_id = c("S1", "S3"),
  expr = c(10.5, 7.8)
)
```

数据只有两列四行的情况下，你很有可能一眼就能发现，样本S2和S4根本没有表达量数据。这可能是因为实验中两个样本的测量出现问题，或者因为某些原因该测量并未进行，而研究者选择了直接清除。当你无视这种情况，你很有可能会在需要整合数据集时出现问题，因此，你至少要提前知道有哪些样本没有出现在真正记录测量值的数据集中。为此你可以这样做：


```R
missing_samples <- df_meta |> dplyr::anti_join(df_expr, by = "sample_id")
missing_samples
```


<table class="dataframe">
<caption>A data.frame: 2 × 2</caption>
<thead>
	<tr><th scope=col>sample_id</th><th scope=col>tissue</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>S2</td><td>Breast</td></tr>
	<tr><td>S4</td><td>Liver </td></tr>
</tbody>
</table>



在样本量巨大的单细胞测序数据集中，这个方法可能会很有用。

## 去除首尾空格的方法:*stringr::str_trim()*

这可能是最有用的函数，你不知道什么时候就会被莫名其妙的空格阴一手。在那些大型数据集中，莫名其妙的空格总是出现在各种奇怪的地方。


```R
df_meta <- data.frame(
  sample_id = c("S1", "S2", "S3"),
  tissue = c("Lung", "Breast", "Colon")
)

df_expr <- data.frame(
  sample_id = c(" S1", "S2 ", "S3"), 
  expr = c(10.5, 8.2, 7.8)
)
```


```R
df_meta
```


<table class="dataframe">
<caption>A data.frame: 3 × 2</caption>
<thead>
	<tr><th scope=col>sample_id</th><th scope=col>tissue</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>S1</td><td>Lung  </td></tr>
	<tr><td>S2</td><td>Breast</td></tr>
	<tr><td>S3</td><td>Colon </td></tr>
</tbody>
</table>




```R
df_expr
```


<table class="dataframe">
<caption>A data.frame: 3 × 2</caption>
<thead>
	<tr><th scope=col>sample_id</th><th scope=col>expr</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> S1</td><td>10.5</td></tr>
	<tr><td>S2 </td><td> 8.2</td></tr>
	<tr><td>S3 </td><td> 7.8</td></tr>
</tbody>
</table>



你说说，看得到S1前面以及S2后面的空格吗？你从csv文件中将数据集导入为数据框，这种问题根本没法用你那俩eyeball Mk.1发现。当你试图用join来整合，你就会发现匹配莫名失败


```R
df_meta |> left_join(df_expr, by = "sample_id") # 这里可不跟你搞模糊数学，匹配不上就是匹配不上，看着一样也不行。
```


<table class="dataframe">
<caption>A data.frame: 3 × 3</caption>
<thead>
	<tr><th scope=col>sample_id</th><th scope=col>tissue</th><th scope=col>expr</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>S1</td><td>Lung  </td><td> NA</td></tr>
	<tr><td>S2</td><td>Breast</td><td> NA</td></tr>
	<tr><td>S3</td><td>Colon </td><td>7.8</td></tr>
</tbody>
</table>



为此，只需要用str_trim()函数清理空格，就可以成功匹配了。


```R
df_expr <- df_expr |>
  mutate(sample_id = stringr::str_trim(sample_id))

df_meta |> left_join(df_expr, by = "sample_id")
```


<table class="dataframe">
<caption>A data.frame: 3 × 3</caption>
<thead>
	<tr><th scope=col>sample_id</th><th scope=col>tissue</th><th scope=col>expr</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>S1</td><td>Lung  </td><td>10.5</td></tr>
	<tr><td>S2</td><td>Breast</td><td> 8.2</td></tr>
	<tr><td>S3</td><td>Colon </td><td> 7.8</td></tr>
</tbody>
</table>



## 把复合编码拆成结构化列:*tidyr::separate()*

我猜你可能见过这种格式的样本列：


```R
df_raw <- data.frame(
  code = c("Lung_001", "Breast_014", "Colon_203")
)
df_raw
```


<table class="dataframe">
<caption>A data.frame: 3 × 1</caption>
<thead>
	<tr><th scope=col>code</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Lung_001  </td></tr>
	<tr><td>Breast_014</td></tr>
	<tr><td>Colon_203 </td></tr>
</tbody>
</table>



这种格式当然很明确地标注了样本的来源，但如果想要根据来源部位来分组，进行差异表达分析这类操作，这种结构就比较麻烦了（当然这种结构可以用正则表达式来解决，但那又是另一个话题了）因此，我们可以这样拆分编码：


```R
df_clean <- df_raw |>
  tidyr::separate(
    col = code,
    into = c("tissue", "sample_id"),
    sep = "_"
  )

df_clean
```


<table class="dataframe">
<caption>A data.frame: 3 × 2</caption>
<thead>
	<tr><th scope=col>tissue</th><th scope=col>sample_id</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Lung  </td><td>001</td></tr>
	<tr><td>Breast</td><td>014</td></tr>
	<tr><td>Colon </td><td>203</td></tr>
</tbody>
</table>



就这样，我们就可以用tissue列来进行分组了。相信我，不是所有的数据集都会给你提供tissue列的，尤其是数据集原作者也没有想到用这种方法进行分组的情况下。

*备注：‘tidyr：：’这种结构就是显式引用，当你引用的多个库有同名函数时最好使用这种引用方法，不然也有可能出现莫名奇妙的问题。当程序出现莫名其妙的bug时，要想到这种可能性。*
