---
layout: post
title: "Tech Log Clean"
math: true
---

# 2026年3月19日 研究日志
今日主题：整理最近学习的生物信息学相关概念、算法及python函数。同时将一些小技巧也记录下来。

### 汉明距离 Hamming distance
在等长字符串之间，对应位置上字符不同的数量。 \
$d_H(p, q) = \sum_{i=1}^{k} [p_i \ne q_i]$ \
在生物信息学中，是比较DNA/RNA/多肽序列相似度最基础的方法,很多问题都需要应用该距离作为解决问题的一部分。


```python
def hamming_distance(p, q):
    """
    Compute the Hamming distance between two strings.

    Parameters
    ----------
    p : str
        The first string to compare.
    q : str
        The second string to compare.

    Returns
    -------
    int
        The number of positions at which the corresponding characters differ.
        Only positions up to the length of the shorter string are compared.
    """·
    return sum(1 for a, b in zip(p, q) if a != b)
```

### 字符串模式计数（Pattern Counting）
统计某个模式串（Pattern）在整个文本串（Text）中出现的次数。在生物信息学中用于统计序列中某个k-mer出现的频率。是进行motif搜索的基础步骤

*k-mer(k-长度子串) 指的是从一条生物序列（DNA/RNA/蛋白质）中截取的长度为𝑘的连续子串。*

*Motif(基序) 在多条DNA或蛋白质序列中反复出现的、具有某种生物学功能的短序列模式，它通常有轻微变异（一直用有突变没办法啦），但核心结构保持一致（不一致导致功能出问题的话个体也很难存活并传递基因）*


```python
def PatternCount(Text, Pattern):
    """
    Count how many times a given pattern appears in a text.

    Parameters
    ----------
    Text : str
        The full string in which the pattern will be searched.
    Pattern : str
        The substring (pattern) to count within the text.

    Returns
    -------
    int
        The total number of occurrences of Pattern in Text,
        including overlapping occurrences.
    """
    count = 0
    pattern_len = len(Pattern)
    
    # Iterate over all possible starting positions in Text
    for i in range(len(Text) - pattern_len + 1):
        # Check whether the substring starting at position i matches Pattern
        if Text[i:i+pattern_len] == Pattern:
            count += 1
    
    return count
```

### 中位串问题 Median String Problem
在字符串集合中寻找一个“中位串”，使得该串到集合中所有其他串的‌总距离最小。
$$
d(p, DNA) = \sum_{s \in DNA} \min_{i} d_H(p, s[i:i+k])
$$
$p$ 候选的k-mer(长度为k的字符串) \
$DNA$ 序列集合，例如 {𝑠1,𝑠2,…,𝑠𝑛} \
对每条序列$𝑠$，我们寻找与$𝑝$最相似的那一段（即最小汉明距离）。然后把所有序列的最小距离加起来。得到的$d(p, DNA)$就是对整个数据集的“代表性”评分。距离越小，说明𝑝越能代表所有序列，也越可能是潜在的 motif。在生物信息学中，对于DNA序列，其调控元件（如转录因子结合位点）往往相对较短且不完全一致（复制多了总要突变的）。因为调控元件往往分布在多条序列中（毕竟调控元件的很多功能大多数序列都要调用），所以找到一个与所有序列都最相似的k-mer，很有可能就是潜在的motif（验证工作交给实验室里的生物学家咯）。

以下是中位串问题的暴力求解法(Brute-force exact solution，就是遍历，只要我把所有可能性都验证一遍，什么东西都能找到，虽然复杂度高，但这是该问题的理论基准)：

↓ d($pattern$, $text$) 计算 $pattern$ 与 $text$ 中所有长度为 $k$ 的子串之间的最小汉明距离


```python
def hamming_distance(p, q):
    return sum(1 for a, b in zip(p, q) if a != b)

def d(pattern, text):
    """
    Compute the minimum Hamming distance between a given pattern and
    any k-length substring of a longer text sequence.

    Parameters
    ----------
    pattern : str
        The k-mer whose similarity to the text will be evaluated.
        Must be a string of length k.
    text : str
        The longer sequence in which all possible k-length windows
        will be compared against the pattern.

    Returns
    -------
    int
        The smallest Hamming distance between `pattern` and any
        substring of `text` with the same length. This value represents
        how closely the pattern matches the best-aligned region of the text.

    Notes
    -----
    This function is a core component of the Median String algorithm.
    It measures the distance between a candidate k-mer and a single DNA
    sequence by finding the most similar window within that sequence.
    """
    k = len(pattern)
    return min(
        hamming_distance(pattern, text[i:i+k])
        for i in range(len(text) - k + 1)
    )
```

↓d_pattern_dna($pattern$, $dna$) 计算$pattern$与所有$DNA$序列的总距离（即对每条序列调用d($pattern$, $text$)）


```python
def d_pattern_dna(pattern, dna):
    """
    Compute the total distance between a given pattern and a collection
    of DNA sequences. The distance to each sequence is defined as the
    minimum Hamming distance between the pattern and any k-length
    substring within that sequence (as computed by d(pattern, text)).

    Parameters
    ----------
    pattern : str
        The candidate k-mer whose overall similarity to the DNA dataset
        will be evaluated.
    dna : list of str
        A list of DNA sequences. Each sequence will be compared against
        the pattern to determine its minimum-distance alignment.

    Returns
    -------
    int
        The sum of distances between `pattern` and each sequence in `dna`.
        This value represents how well the pattern matches the entire
        collection of sequences. Lower values indicate a better match.

    Notes
    -----
    This function is used in the Median String algorithm to evaluate
    how good a candidate k-mer is across all sequences in the dataset.
    """
    return sum(d(pattern, seq) for seq in dna)
```

↓核心算法：枚举所有可能的 k-mer，计算它们的总距离，返回距离最小的那个。


```python
def median_string(dna, k):
    """
    Find the k-mer (pattern of length k) that minimizes the total distance
    to a collection of DNA sequences. This is the classical Median String
    Problem in bioinformatics.

    The function exhaustively enumerates all possible k-mers over the
    alphabet {A, C, G, T}, computes their total distance to the DNA
    dataset using d_pattern_dna(), and returns the k-mer with the
    smallest total distance.

    Parameters
    ----------
    dna : list of str
        A list of DNA sequences over which the median string will be
        computed.
    k : int
        The length of the k-mer to search for.

    Returns
    -------
    str
        The k-mer that minimizes the total distance to all sequences
        in `dna`. If multiple k-mers achieve the same minimum distance,
        the first encountered in lexicographic order is returned.

    Notes
    -----
    This brute-force implementation has time complexity O(4^k * n * k),
    where n is the number of sequences. It is exact but computationally
    expensive for large k.
    """
    best_pattern = None
    best_distance = float("inf")

    for pattern_tuple in product("ACGT", repeat=k):
        pattern = "".join(pattern_tuple)
        distance = d_pattern_dna(pattern, dna)

        if distance < best_distance:
            best_distance = distance
            best_pattern = pattern

    return best_pattern
```

### 编写python代码的小技巧

↓小技巧：使用上下文管理语句with打开txt格式序列文件。比手动打开和关闭更安全更简洁。


```python
text_file = 'path/to/file'
with open('text_file','r', encoding = 'utf-8') as f:
    content = f.read()
```

↓小技巧：原始字符串前缀r"" 反斜杠不再被当作转义符，符合windows路径习惯。而且用到正则表达式的时候也是必不可少（说来惭愧，我自己正则表达式也用的稀里糊涂）


```python
path = r"C:\new\folder"
```
