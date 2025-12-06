# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.gridspec as gridspec
import seaborn as sns
from Bio import SeqIO
import os
from collections import defaultdict
import matplotlib as mpl
from matplotlib.font_manager import FontProperties
import math

# %%
# 设置全局样式
plt.rcParams.update({
    'font.family': 'Arial',
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
})

# 1. 解析输入文件
fasta_file = '../results/filtered_ir_proteins.fasta'
pfam_file = '../results/pfamscan_results_new_3domain.txt'

def parse_pfam_results(pfam_file, target_domains):
    """
    解析Pfam扫描结果文件
    返回一个字典：{序列ID: 结构域信息列表}
    """
    domain_df = pd.read_csv(pfam_file, sep='\s+', comment='#', header=None,
                           names=['seq_id', 'aln_start', 'aln_end', 'env_start', 'env_end', 
                                  'hmm_acc', 'hmm_name', 'type', 'hmm_start', 'hmm_end', 
                                  'hmm_len', 'bit_score', 'evalue', 'significance', 'clan'])
    domain_df = domain_df[domain_df['hmm_name'].isin(target_domains)]
    grouped = domain_df.groupby('seq_id')
    domain_dict = {}
    for seq_id, group in grouped:
        domains = []
        for _, row in group.iterrows():
            domains.append({
                'type': row['hmm_name'],
                'start': row['env_start'],
                'end': row['env_end'],
                'bit_score': row['bit_score'],
                'evalue': row['evalue'],
                'length': row['env_end'] - row['env_start'] + 1
            })
        domain_dict[seq_id] = sorted(domains, key=lambda x: x['start'])
    return domain_dict

# 提取PAS和GAF的前30条
pas_gaf_domains = ['PAS', 'GAF']
pas_gaf_full = parse_pfam_results(pfam_file, pas_gaf_domains)
pas_gaf_dict = {k: pas_gaf_full[k] for k in list(pas_gaf_full.keys())[:30]} # 

# 提取PHY的全部
phy_domains = ['PHY']
phy_dict = parse_pfam_results(pfam_file, phy_domains)

# 合并到一个domain_dict
domain_dict = {}
domain_dict.update(pas_gaf_dict)
domain_dict.update(phy_dict)

def get_sequence_lengths(fasta_file):
    """获取每个序列的长度"""
    seq_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_lengths[record.id] = len(record.seq)
    return seq_lengths

# 获取序列长度
seq_lengths = get_sequence_lengths(fasta_file)

print(f"分析完成: {len(domain_dict)}条序列包含PAS/GAF/PHY结构域")

# %%
# 新增：绘制PAS/GAF/PHY数量柱状图
def plot_domain_counts(domain_dict, output_path='../results/domain_counts.pdf', figsize=(4, 6)):
    domain_colors = {'PAS': '#FF6B6B', 'GAF': '#4ECDC4', 'PHY': '#FFD166'}
    # 统计各类型结构域数量
    counts = {'PAS': 0, 'GAF': 0, 'PHY': 0}
    for domains in domain_dict.values():
        for d in domains:
            if d['type'] in counts:
                counts[d['type']] += 1
    # 绘图
    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.bar(
        list(counts.keys()), 
        list(counts.values()), 
        color=[domain_colors[k] for k in counts.keys()], 
        width=0.8
    )
    # 设置字体大小
    ax.set_ylabel('Count', fontsize=16)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    # 去掉x轴标题和总标题
    # ax.set_xlabel('Domain Type', fontsize=16)
    # ax.set_title('Number of PAS, GAF, PHY Domains', fontsize=18)
    ax.bar_label(bars, fontsize=16)
    # 去除上、右边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # 只保留xy轴
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    plt.tight_layout()
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    return fig, ax

# 调用绘图并保存为PDF
plot_domain_counts(domain_dict)

# # %%
# def plot_domain_architectures(domain_dict, seq_lengths, figsize=(8, 8)):
#     """
#     绘制结构域架构图：显示每个序列中结构域的位置和顺序
#     """
#     # 定义结构域颜色映射
#     domain_colors = {
#         'PAS': '#FF6B6B',  # 红色系
#         'GAF': '#4ECDC4',  # 青色系
#         'PHY': '#FFD166'   # 黄色系
#     }
    
#     # 创建图形和坐标轴
#     fig = plt.figure(figsize=figsize)
#     gs = gridspec.GridSpec(2, 1, height_ratios=[1.5, 0.1])
#     ax1 = plt.subplot(gs[0])  # 主图
#     ax2 = plt.subplot(gs[1])  # 图例
    
#     # 准备数据：按序列ID排序
#     seq_ids = sorted(domain_dict.keys())
    
#     # 计算每个序列的最大X范围
#     max_length = max(seq_lengths.values())
#     max_length = 900
    
#     # 为每个序列绘制结构域
#     patch_list = []
#     y_pos = 0
    
#     # 添加每个序列的域块
#     for i, seq_id in enumerate(seq_ids):
#         domains = domain_dict[seq_id]
#         seq_length = seq_lengths[seq_id]
        
#         # 绘制序列背景
#         ax1.add_patch(Rectangle((0, y_pos-0.2), seq_length, 0.8, 
#                               facecolor='#F0F0F0', edgecolor='gray', alpha=0.5))
        
#         # 绘制每个结构域
#         for domain in domains:
#             start = domain['start']
#             end = domain['end']
#             domain_type = domain['type']
#             domain_length = end - start + 1
            
#             patch = Rectangle((start, y_pos-0.1), domain_length, 0.6,
#                              facecolor=domain_colors[domain_type], edgecolor='black', alpha=0.8)
#             ax1.add_patch(patch)
            
#             # 在域中间添加文本标签（仅当域长度足够时），向上偏移0.3
#             if domain_length > 20:
#                 ax1.text((start+end)/2, y_pos+0.2, domain_type, 
#                         ha='center', va='center', fontsize=9, fontweight='bold')
                
#             patch_list.append(patch)
        
#         # 添加序列ID标签（右对齐），向上偏移0.3
#         ax1.text(seq_length + max_length*0.01, y_pos+0.3, seq_id, 
#                 ha='left', va='center', fontsize=8, color='black')
        
#         # 添加序列长度刻度，向下偏移0.2
#         ax1.text(seq_length, y_pos-0.2, f"{seq_length} aa", 
#                 ha='right', va='top', fontsize=7, color='gray')
        
#         y_pos -= 1  # 下移一行
    
#     # 设置Y轴：每个序列一行
#     yticks = np.arange(-(len(seq_ids)-1), 1, 1)[::-1]
#     ax1.set_yticks(yticks)
#     # 将y轴标签向上偏移0.1
#     for tick, label in zip(yticks, seq_ids):
#         ax1.text(-max_length*0.03, tick+0.2, label, ha='right', va='center', fontsize=10, color='black', clip_on=False)
#     ax1.set_yticklabels([])  # 不显示默认标签
    
#     # 设置X轴
#     ax1.set_xlim(0, max_length*1.1)
#     ax1.set_xlabel('Amino Acid Position')
    
#     # 添加标题
#     ax1.set_title('Domain Architecture of Infrared-Optogenetic Proteins', fontsize=16, pad=20)
    
#     # 添加网格线
#     ax1.grid(True, which='major', axis='x', linestyle='--', color='gray', alpha=0.3)
    
#     # 创建图例
#     legend_patches = []
#     for domain_type, color in domain_colors.items():
#         legend_patches.append(
#             Rectangle((0,0), 1, 1, facecolor=color, edgecolor='black', alpha=0.8, label=domain_type)
#         )
    
#     # 在底部添加图例
#     ax2.axis('off')
#     legend = ax2.legend(handles=legend_patches, loc='center', 
#                         ncol=2, frameon=False, fontsize=14)
#     legend.set_title('Domain Types', prop={'weight': 'bold', 'size': 14})
    
#     plt.tight_layout()
#     output_plot = '../results/infrared_optogenetic_domain_architectures'
#     plt.savefig(f"{output_plot}.pdf", format="pdf", bbox_inches="tight")
#     return fig, ax1

# # 2. 创建分析报告
# report_fig = plt.figure(figsize=(8, 8))

# # 结构域架构图
# domain_arch_fig, _ = plot_domain_architectures(domain_dict, seq_lengths, figsize=(8, 6))

# # %%

# # %%
# def plot_domain_density(domain_dict, seq_lengths, figsize=(5, 8)):
#     """
#     绘制结构域分布密度图：显示PAS/GAF/PHY在序列中出现的频次分布
#     """
#     # 创建位置数据列表
#     density_data = {'Position': [], 'Domain': []}
    
#     max_length = max(seq_lengths.values())
#     bin_size = 50  # 每50个氨基酸为一个区间

#     # 收集所有结构域位置
#     for seq_id, domains in domain_dict.items():
#         for domain in domains:
#             # 添加每个位置点
#             density_data['Position'].append(domain['start'])
#             density_data['Domain'].append(domain['type'])
#             density_data['Position'].append(domain['end'])
#             density_data['Domain'].append(domain['type'])
    
#     # 转换为DataFrame
#     df = pd.DataFrame(density_data)
    
#     # 创建分面网格图
#     fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True)
#     # domain_colors = {'PAS': '#FF6B6B', 'GAF': '#4ECDC4', 'PHY': '#FFD166'}
#     domain_colors = {'GAF': '#4ECDC4', 'PHY': '#FFD166'}

#     # for i, domain_type in enumerate(['PAS', 'GAF', 'PHY']): 
#     for i, domain_type in enumerate(['GAF', 'PHY']):
#         ax = axes[i]
        
#         # 提取特定结构域的数据
#         subset = df[df['Domain'] == domain_type]
        
#         # 绘制KDE密度曲线
#         sns.kdeplot(data=subset, x='Position', fill=True, alpha=0.5, 
#                    color=domain_colors[domain_type], ax=ax, label=f'{domain_type} density')
        
#         # 绘制直方图
#         sns.histplot(data=subset, x='Position', bins=range(0, 1000+bin_size, bin_size),
#                     stat='percent', alpha=0.7, color=domain_colors[domain_type], 
#                     element='step', ax=ax)
        
#         # 设置中位线
#         median_pos = np.median(subset['Position'])
#         ax.axvline(median_pos, color='black', linestyle='--', linewidth=1)
#         ax.text(median_pos+10, ax.get_ylim()[1]*0.9, 
#                f'Median: {int(median_pos)}', fontsize=9)
        
#         ax.set_title(f'{domain_type} Domain Position Distribution')
#         ax.set_ylabel('Frequency (%)')
#         ax.set_xlim(0, 1000)  # 缩小x轴范围为0-1000
    
#     # 设置公共X轴标签
#     axes[-1].set_xlabel('Amino Acid Position')
#     plt.suptitle('Positional Distribution of Key Domains', fontsize=16, y=1)
#     plt.tight_layout()
#     return fig
# # %%
# def plot_domain_statistics(domain_dict, figsize=(14, 10)):
#     """
#     绘制结构域统计信息：大小、位点得分、E-value分布
#     """
#     # 准备数据
#     stats_data = {'Domain': [], 'BitScore': [], 'Evalue': [], 'Length': []}
    
#     for domains in domain_dict.values():
#         for domain in domains:
#             stats_data['Domain'].append(domain['type'])
#             stats_data['BitScore'].append(domain['bit_score'])
#             stats_data['Evalue'].append(domain['evalue'])
#             stats_data['Length'].append(domain['length'])
    
#     stats_df = pd.DataFrame(stats_data)
    
#     # 初始化画布
#     fig = plt.figure(figsize=figsize)
#     gs = gridspec.GridSpec(2, 2)
    
#     # 1. 结构域长度分布 (箱线图和小提琴图)
#     ax1 = plt.subplot(gs[0, 0])
#     sns.boxplot(data=stats_df, x='Domain', y='Length', palette=['#FF6B6B', '#4ECDC4', '#FFD166'], ax=ax1)
#     sns.stripplot(data=stats_df, x='Domain', y='Length', color='black', alpha=0.4, ax=ax1, jitter=True)
#     ax1.set_title('Domain Length Distribution')
#     ax1.set_xlabel('')
#     ax1.set_ylabel('Amino Acids')
    
#     # 2. Bit Score分布 (小提琴图)
#     ax2 = plt.subplot(gs[0, 1])
#     sns.violinplot(data=stats_df, x='Domain', y='BitScore', 
#                   inner='quartile', palette=['#FF6B6B', '#4ECDC4', '#FFD166'], ax=ax2)
#     ax2.set_title('Bit Score Distribution (Domain Quality)')
#     ax2.set_ylabel('Bit Score')
    
#     # 3. E-value分布 (对数分布)
#     ax3 = plt.subplot(gs[1, 0])
#     for domain in ['PAS', 'GAF', 'PHY']:
#         subset = stats_df[stats_df['Domain'] == domain]
#         sns.histplot(np.log10(subset['Evalue'] + 1e-100), kde=True, 
#                     label=domain, alpha=0.7, ax=ax3,
#                     color={'PAS': '#FF6B6B', 'GAF': '#4ECDC4', 'PHY': '#FFD166'}[domain])
    
#     ax3.set_title('E-value Distribution (log scale)')
#     ax3.set_xlabel('log$_{10}$(E-value)')
#     ax3.set_ylabel('Count')
#     ax3.legend(title='Domain')
    
#     # 4. 结构域数量关系
#     ax4 = plt.subplot(gs[1, 1])
#     count_df = stats_df.groupby(['Domain']).size().reset_index(name='Count')
#     count_df['Percent'] = count_df['Count'] / count_df['Count'].sum() * 100
    
#     # 环形图
#     wedges, texts, autotexts = ax4.pie(
#         count_df['Count'], 
#         labels=count_df['Domain'],
#         colors=['#FF6B6B', '#4ECDC4', '#FFD166'],
#         autopct='%1.1f%%',
#         startangle=90,
#         wedgeprops={'edgecolor': 'w', 'linewidth': 2},
#         textprops={'fontsize': 12, 'fontweight': 'bold'},
#         pctdistance=0.85
#     )
    
#     # 添加标题
#     circle = plt.Circle((0,0), 0.7, fc='white')
#     ax4.add_artist(circle)
#     ax4.set_aspect('equal')
#     ax4.set_title('Domain Type Distribution', fontsize=14)
    
#     plt.suptitle('Domain Statistics for Infrared-Optogenetic Proteins', fontsize=16, y=0.95)
#     plt.tight_layout(pad=3.0)
#     return fig
# # %%
# def plot_domain_combination_network(domain_dict, figsize=(12, 12)):
#     """
#     绘制结构域组合网络图：显示结构域在序列中的组合关系
#     """
#     # 生成结构域组合模式
#     patterns = {}
#     for seq_id, domains in domain_dict.items():
#         # 创建结构域序列表示（按位置排序）
#         domain_sequence = '-'.join(sorted(d['type'] for d in domains))
#         patterns[seq_id] = domain_sequence
    
#     # 统计每种组合的频率
#     from collections import Counter
#     pattern_count = Counter(patterns.values())
    
#     # 创建图结构
#     import networkx as nx
#     G = nx.Graph()
    
#     # 添加节点（每种结构域）
#     domain_types = ['PAS', 'GAF', 'PHY']
#     for domain in domain_types:
#         G.add_node(domain, node_type='domain')
    
#     # 添加组合模式节点
#     for pattern, count in pattern_count.items():
#         # 只包含所有三种结构域的组合
#         if all(d in pattern for d in domain_types):
#             G.add_node(pattern, node_type='pattern', size=count*10)
    
#     # 添加边：结构域节点到组合节点
#     for pattern in pattern_count:
#         # 只包含所有三种结构域的组合
#         if all(d in pattern for d in domain_types):
#             for domain in pattern.split('-'):
#                 if domain in domain_types and pattern in G.nodes:
#                     G.add_edge(domain, pattern, weight=pattern_count[pattern])
    
#     # 创建位置布局
#     pos = {}
    
#     # 结构域节点在中心周围等距排列
#     angle_step = 2 * math.pi / len(domain_types)
#     radius = 3
#     for i, domain in enumerate(domain_types):
#         angle = i * angle_step
#         pos[domain] = (radius * math.cos(angle), radius * math.sin(angle))
    
#     # 组合节点在中心
#     for pattern in pattern_count:
#         if pattern in G.nodes and G.nodes[pattern]['node_type'] == 'pattern':
#             pos[pattern] = (0, 0)
    
#     # 绘制网络图
#     fig, ax = plt.subplots(figsize=figsize)
    
#     # 定义节点颜色
#     domain_colors = {'PAS': '#FF6B6B', 'GAF': '#4ECDC4', 'PHY': '#FFD166'}
#     node_colors = []
    
#     # 绘制节点
#     for node in G.nodes():
#         if G.nodes[node]['node_type'] == 'domain':
#             node_colors.append(domain_colors[node])
#         else:
#             node_colors.append('purple')
    
#     # 绘制边
#     nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.4, ax=ax)
    
#     # 绘制节点
#     nx.draw_networkx_nodes(
#         G, pos, node_size=[G.nodes[n].get('size', 800) for n in G.nodes()],
#         node_color=node_colors, alpha=0.8, edgecolors='black', ax=ax
#     )
    
#     # 添加标签
#     node_labels = {}
#     for node in G.nodes():
#         if G.nodes[node]['node_type'] == 'domain':
#             node_labels[node] = node
#         else:
#             # 显示组合模式和序列数量
#             node_labels[node] = f"{node}\n({pattern_count[node]} seqs)"
    
#     nx.draw_networkx_labels(
#         G, pos, labels=node_labels, 
#         font_size=10, font_weight='bold', ax=ax
#     )
    
#     # 添加标题
#     plt.title('Domain Combination Network', fontsize=16)
#     plt.axis('off')
    
#     return fig

# # 主分析流程
# if __name__ == "__main__":
#     # 1. 解析输入文件
#     fasta_file = '../data/results/filtered_ir_proteins.fasta'
#     pfam_file = '../data/results/pfamscan_results_new_3domain.txt'
    
#     # 解析Pfam结果
#     domain_dict = parse_pfam_results(pfam_file)
    
#     # 获取序列长度
#     seq_lengths = get_sequence_lengths(fasta_file)
    
#     print(f"分析完成: {len(domain_dict)}条序列包含PAS/GAF/PHY结构域")
    
#     # 2. 创建分析报告
#     report_fig = plt.figure(figsize=(18, 20))
    
#     # 结构域架构图
#     domain_arch_fig, _ = plot_domain_architectures(domain_dict, seq_lengths, figsize=(18, 12))
    
#     # 结构域密度图
#     density_fig = plot_domain_density(domain_dict, seq_lengths, figsize=(16, 8))
    
#     # 结构域统计图
#     stats_fig = plot_domain_statistics(domain_dict, figsize=(16, 10))
    
#     # 结构域组合网络
#     network_fig = plot_domain_combination_network(domain_dict, figsize=(12, 12))
    
#     # 3. 保存所有图表
#     domain_arch_fig.savefig('infrared_optogenetic_domain_architectures.png', dpi=300, bbox_inches='tight')
#     density_fig.savefig('domain_position_density.png', dpi=300, bbox_inches='tight')
#     stats_fig.savefig('domain_statistics.png', dpi=300, bbox_inches='tight')
#     network_fig.savefig('domain_combination_network.png', dpi=300, bbox_inches='tight')
    
#     # 4. 可选：创建PDF报告
#     print("可视化完成! 图表已保存到当前目录")
# %%
def plot_domain_architectures(domain_dict, seq_lengths, figsize):
    """
    绘制结构域架构图：显示每个序列中结构域的位置和顺序
    """
    # 定义结构域颜色映射
    domain_colors = {
        'PAS': '#FF6B6B',  # 红色系
        'GAF': '#4ECDC4',  # 青色系
        'PHY': '#FFD166'   # 黄色系
    }
    
    # 创建图形和坐标轴
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 1)
    ax1 = plt.subplot(gs[0])  # 主图
    
    # 准备数据：按序列ID排序
    seq_ids = sorted(domain_dict.keys())
    
    # 计算每个序列的最大X范围
    max_length = max(seq_lengths.values())
    max_length = min(max_length, 1400)  # 图的长度减小
    
    # 为每个序列绘制结构域
    patch_list = []
    y_gap = 0.3  # 增大y轴间隔
    y_pos = 0
    
    # 添加每个序列的域块
    for i, seq_id in enumerate(seq_ids):
        domains = domain_dict[seq_id]
        seq_length = seq_lengths[seq_id]
        
        # 绘制序列背景（瘦一点）
        ax1.add_patch(Rectangle((0.05, y_pos-0.05), seq_length, 0.25, 
                              facecolor='#F0F0F0', edgecolor='gray', alpha=0.5))
        
        # 绘制每个结构域（瘦一点）
        for domain in domains:
            start = domain['start']
            end = domain['end']
            domain_type = domain['type']
            domain_length = end - start + 0.1
            
            patch = Rectangle((start, y_pos), domain_length, 0.15,
                             facecolor=domain_colors[domain_type], edgecolor='black', alpha=0.8)
            ax1.add_patch(patch)
            
            # 在域中间添加文本标签（仅当域长度足够时），向上偏移0.18
            if domain_length > 20:
                ax1.text((start+end)/2, y_pos+0.05, domain_type, 
                        ha='center', va='center', fontsize=16, fontweight='bold')
                
            patch_list.append(patch)
        
        # 添加序列长度刻度，向下偏移0.1
        ax1.text(seq_length + max_length*0.01, y_pos+0.05, f"{seq_length} aa", 
                ha='left', va='center', fontsize=16, color='gray')
        
        y_pos -= y_gap  # 下移一行，间隔更大
    
    # 设置Y轴：每个序列一行
    yticks = np.arange(0, -y_gap*len(seq_ids), -y_gap)
    ax1.set_yticks(yticks)
    # 将y轴标签向上偏移0.1
    for tick, label in zip(yticks, seq_ids):
        ax1.text(-max_length*0.03, tick+0.1, label, ha='right', va='center', fontsize=16, color='black', clip_on=False)
    ax1.set_yticklabels([])  # 不显示默认标签
    
    # 设置X轴
    ax1.set_xlim(0, max_length*1.05)
    ax1.set_xlabel('Amino Acid Position', fontsize=16)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    
    # 不添加大标题
    
    # 添加网格线
    ax1.grid(True, which='major', axis='x', linestyle='--', color='gray', alpha=0.3)
    
    # 上下各保留0.5个y_gap的空白
    ax1.set_ylim(yticks[-1] - y_gap*0.5, yticks[0] + y_gap*0.9)
    
    # 创建图例
    legend_patches = []
    for domain_type, color in domain_colors.items():
        legend_patches.append(
            Rectangle((0,0), 1, 1, facecolor=color, edgecolor='black', alpha=0.8, label=domain_type)
        )
    
    # 在主图右下角添加图例
    legend = ax1.legend(handles=legend_patches, loc='lower right', 
                        ncol=1, frameon=False, fontsize=16, title='Domain Types', title_fontsize=16)
    
    plt.tight_layout()
    output_plot = '../results/infrared_optogenetic_domain_architectures'
    plt.savefig(f"{output_plot}.pdf", format="pdf", bbox_inches="tight")
    return fig, ax1

# 2. 创建分析报告
report_fig = plt.figure(figsize=(8, 10))

# 结构域架构图
domain_arch_fig, _ = plot_domain_architectures(domain_dict, seq_lengths, figsize=(12, 20))

# %%
