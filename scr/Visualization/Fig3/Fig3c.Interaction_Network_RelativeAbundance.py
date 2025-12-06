# %%
import pandas as pd
import numpy as np
# from skbio.stats.composition import clr
import networkx as nx
import community  # python-louvain包
from sklearn.ensemble import RandomForestClassifier
import os

# %%
# 创建输出目录
if not os.path.exists("network_data"):
    os.makedirs("network_data")

# 设置工作环境
np.random.seed(123)  # 设置随机种子，保证结果可重复

# --------------------
# 1. 数据导入与预处理
# --------------------

# 设置工作路径
data_dir = '../../data/Fig3/Network/'  # 更改为您的数据目录

# 文件路径
abundance_file = os.path.join(data_dir, 'metagenoNR_table.umap_df.cluster3.csv')
taxonomy_file =  os.path.join(data_dir, 'taxonomy.csv')
metadata_file = os.path.join(data_dir, 'metadata.csv')

# 步骤1：数据读取与预处理
# 读取ASV表（行为ASV，列为样本）
asv_table = pd.read_csv(abundance_file, index_col=0)

# 读取分类信息（第一列是ASV，第二列是分类字符串）
taxonomy = pd.read_csv(taxonomy_file, index_col=0, header=0)
taxonomy.columns = ['Taxonomy']  # 确保列名为Taxonomy
# %%
# 解析分类字符串为各个分类等级
def parse_taxonomy(tax_string):
    ranks = ['Phylum', 'Genus']
    # ranks = ['Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    split_tax = tax_string.split(';')
    parsed = {}
    for i, tax_level in enumerate(split_tax):
        if i < len(ranks):
            # 移除前缀如'd__'
            if '__' in tax_level:
                # 分离前缀和实际名称
                prefix, name = tax_level.split('__', 1) if tax_level.count('__') == 1 else ('', tax_level)
                parsed[ranks[i]] = name
            else:
                parsed[ranks[i]] = tax_level
    return parsed

# 应用解析函数
taxonomy_parsed = taxonomy['Taxonomy'].apply(parse_taxonomy).apply(pd.Series)
taxonomy_parsed.index = taxonomy.index  # 保持索引一致

# 读取元数据
metadata = pd.read_csv(metadata_file, index_col="Sample")  # 索引为样本名称

# 过滤低丰度ASV（在20%样本中出现,此处可选）
min_samples = int(0.1 * asv_table.shape[1])  # 列数为样本数
asv_table_filtered = asv_table.loc[(asv_table > 0).sum(axis=1) >= min_samples]

# 相对丰度标准化（按列，即按样本求和）
rel_asv = asv_table_filtered.div(asv_table_filtered.sum(axis=0)).T  # 转置为样本×ASV

# %%
# 步骤2：构建相关网络（使用CLR转换的Pearson相关）
def clr_correlation(data):
    # 处理零值
    data = data + 1e-10
    clr_data = np.log(data).apply(lambda x: x - np.mean(x), axis=1)
    return clr_data.corr()

# 计算相关性矩阵
corr_matrix = clr_correlation(rel_asv)

# 设置阈值构建邻接矩阵（绝对值大于0.6），保留相关系数作为权重
adj_matrix = corr_matrix.copy()
adj_matrix[np.abs(adj_matrix) <= 0.6] = 0
np.fill_diagonal(adj_matrix.values, 0)  # 对角线为0

# 移除负权重（Louvain算法不支持负权重）
adj_matrix[adj_matrix < 0] = 0

# 创建图，保留权重
G = nx.from_pandas_adjacency(adj_matrix)

# 步骤3：网络拓扑分析
# 计算中心性指标
topology_metrics = pd.DataFrame({
    'degree': nx.degree_centrality(G),
    'betweenness': nx.betweenness_centrality(G),
    'closeness': nx.closeness_centrality(G)
}).reset_index().rename(columns={'index': 'ASV'})

# 特征向量中心性（可能不收敛，需处理）
try:
    eigenvector = nx.eigenvector_centrality(G, max_iter=1000)
    topology_metrics['eigenvector'] = topology_metrics['ASV'].map(eigenvector)
except nx.PowerIterationFailedConvergence:
    print("Eigenvector centrality failed to converge, using degree centrality instead.")
    topology_metrics['eigenvector'] = topology_metrics['degree']
# %%
# 步骤4：基石物种识别
# Louvain算法划分模块
partition = community.best_partition(G, random_state=42)

# 自动寻找最佳相关性阈值r，使网络既有模块结构又有足够的连接
def find_best_r(corr_matrix, r_range=np.arange(0.3, 0.8, 0.05), min_nodes=10, min_edges=20):
    best_r = None
    best_modularity = -1
    best_G = None
    best_partition = None
    for r in r_range:
        adj = corr_matrix.copy()
        adj[np.abs(adj) <= r] = 0
        np.fill_diagonal(adj.values, 0)
        adj[adj < 0] = 0
        G_tmp = nx.from_pandas_adjacency(adj)
        if G_tmp.number_of_nodes() < min_nodes or G_tmp.number_of_edges() < min_edges:
            continue
        part = community.best_partition(G_tmp, random_state=42)
        modularity = community.modularity(part, G_tmp)
        if modularity > best_modularity:
            best_modularity = modularity
            best_r = r
            best_G = G_tmp
            best_partition = part
    return best_r, best_G, best_partition

# 自动选择r
best_r, G, partition = find_best_r(corr_matrix)
print(f"Best correlation threshold r: {best_r}")

# 计算每个节点的Zi和Pi，并自动寻找最佳阈值
def calculate_zi_pi(G, partition):
    zi_pi = {}
    zi_list = []
    pi_list = []
    for node in G.nodes():
        same_module_nodes = [n for n in G.neighbors(node) if partition.get(n) == partition.get(node)]
        other_module_nodes = [n for n in G.neighbors(node) if partition.get(n) != partition.get(node)]
        ki = len(list(G.neighbors(node)))
        kis = len(same_module_nodes)
        kio = len(other_module_nodes)
        same_module_nodes_in = same_module_nodes.copy()
        same_module_nodes_in.append(node)
        if same_module_nodes_in:
            subgraph = G.subgraph(same_module_nodes_in)
            degrees = [d for _, d in subgraph.degree()]
            k_s_i = np.mean(degrees)
            sigma_k_s = np.std(degrees) if len(degrees) > 1 else 0
        else:
            k_s_i = 0
            sigma_k_s = 0
        if sigma_k_s != 0:
            zi = (kis - k_s_i) / sigma_k_s
        else:
            zi = 0
        if ki > 0:
            pi = kio / ki
        else:
            pi = 0
        zi_list.append(zi)
        pi_list.append(pi)
        zi_pi[node] = {'Zi': zi, 'Pi': pi}
    # 自动用分位数确定阈值
    zi_cut = np.percentile(zi_list, 90)  # v2: 上四分位(75)；v1: 90分位 
    pi_cut = np.percentile(pi_list, 90)
    for node in zi_pi:
        zi = zi_pi[node]['Zi']
        pi = zi_pi[node]['Pi']
        role = "Peripheral"
        if zi > zi_cut and pi > pi_cut:
            role = "Network Hub"
        elif zi > zi_cut:
            role = "Module Hub"
        elif pi > pi_cut:
            role = "Connector"
        zi_pi[node]['Role'] = role
    print(f"Zi threshold: {zi_cut:.2f}, Pi threshold: {pi_cut:.2f}")
    return zi_pi


zi_pi_df = pd.DataFrame(calculate_zi_pi(G, partition)).T.reset_index().rename(columns={'index': 'ASV'})

# %%
# 随机森林验证（使用分组信息）
if 'Group' in metadata.columns:
    rf = RandomForestClassifier(n_estimators=500, random_state=42)
    rf.fit(rel_asv, metadata.loc[rel_asv.index, 'Group'])  # 确保元数据顺序与rel_asv一致
    importances = pd.DataFrame({
        'ASV': rel_asv.columns,
        'RF_Importance': rf.feature_importances_
    })
else:
    importances = pd.DataFrame({
        'ASV': rel_asv.columns,
        'RF_Importance': np.nan
    })
    print("No 'Group' column found in metadata. Skipping randomForest analysis.")

# 合并结果
keystone_species = (
    topology_metrics
    .merge(zi_pi_df, on='ASV', how='left')
    .merge(importances, on='ASV', how='left')
    .merge(taxonomy_parsed, left_on='ASV', right_index=True, how='left')
)

# 识别基石物种
keystone_species['is_keystone'] = False

# 确保所有列都存在
if 'betweenness' in keystone_species and 'Role' in keystone_species and 'RF_Importance' in keystone_species:
    keystone_condition = (
        keystone_species['Role'].isin(['Network Hub', 'Module Hub', 'Connector']) &
        (keystone_species['betweenness'] > keystone_species['betweenness'].quantile(0.6))
    )
    
    # 如果有RF重要性则加入条件
    if not keystone_species['RF_Importance'].isnull().all():
        keystone_condition &= (keystone_species['RF_Importance'] > keystone_species['RF_Importance'].quantile(0.6))
    
    keystone_species['is_keystone'] = keystone_condition

# 保存边列表（包含权重，保留正负相关）
# 重新从原始相关矩阵（corr_matrix）生成边列表，保留正负权重
edges = []
for i in corr_matrix.index:
    for j in corr_matrix.columns:
        if i < j and np.abs(corr_matrix.loc[i, j]) > 0.6:
            edges.append({
                'source': i,
                'target': j,
                'weight': corr_matrix.loc[i, j]
            })
edge_list = pd.DataFrame(edges)
edge_list.to_csv("network_data/network_edges.csv", index=False)

# 保存节点属性
nodes_with_attributes = keystone_species.copy()
nodes_with_attributes['module'] = nodes_with_attributes['ASV'].map(partition)
nodes_with_attributes.to_csv("network_data/network_nodes.csv", index=False)

# 保存基石物种列表
keystone_species[keystone_species['is_keystone']].to_csv("network_data/keystone_species.csv", index=False)

# 保存元数据
metadata.to_csv("network_data/metadata.csv")

print("Python analysis completed. Data saved in 'network_data' directory.")

