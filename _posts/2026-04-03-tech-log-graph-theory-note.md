# 2026年4月3日 研究日志
今日研究内容：学习有关图论和深度优先遍历的基础知识。

说真的，基因可以被理解为图我是今天才知道的，序列分析可以通过将其视作顶点为k-mer，边则是重叠关系的De Brujin图，单细胞组学的kNN也是将细胞视作顶点，相似度视作边。以前一直不求甚解地用着这些做项目，没想到我的工作底层竟然在这里，现在学习了一点，感觉眼界又开阔了一些。

## DFS（深度优先搜索）
从一个节点开始，沿着一条路径尽可能深入 → 找到所有可能的路径


```python
import networkx as nx

def dfs(graph, start):
    visited = set()
    order = []

    def _dfs(node):
        visited.add(node)
        order.append(node)
        for neighbor in graph.neighbors(node):
            if neighbor not in visited:
                _dfs(neighbor)

    _dfs(start)
    return order

G = nx.Graph()
G.add_edges_from([
    ("A", "B"), ("A", "D"),
    ("B", "C"),
    ("C", "E"),
    ("D", "E")
])

print("DFS 遍历顺序:", dfs(G, "A"))
```

    DFS 遍历顺序: ['A', 'B', 'C', 'E', 'D']
    

在基因调控网络（GRN）或 PPI 网络中，DFS 常用于寻找调控链、探索深层路径、检测反馈环（feedback loops），以及在单细胞轨迹树中遍历分化路径。

## BFS（广度优先搜索）
按层扩散，用于找到最短路径或局部结构


```python
import networkx as nx
from collections import deque

def bfs(graph, start):
    visited = set([start])
    queue = deque([start])
    order = []

    while queue:
        node = queue.popleft()
        order.append(node)
        for neighbor in graph.neighbors(node):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)

    return order

G = nx.Graph()
G.add_edges_from([
    ("A", "B"), ("A", "D"),
    ("B", "C"),
    ("C", "E"),
    ("D", "E")
])

print("BFS 遍历顺序:", bfs(G, "A"))
```

    BFS 遍历顺序: ['A', 'B', 'D', 'C', 'E']
    

BFS 在组学中常用于寻找最短生物通路（如代谢网络中的最短反应链）、单细胞 kNN 图中的连通性检测、网络扩散分析（disease gene propagation），以及识别距离某个基因最近的功能模块。

*说来惭愧，我现在还对这些算法懵懵懂懂，所以没法教学，今天就做个展示，之后再深入了解吧。*
