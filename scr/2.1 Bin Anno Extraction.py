# === 
# This scripts will extract certain annotation (e.g. CO2 fixation related KEGG level 3 pathways) from the annotation file
# ===
# """
# Input
# - Anno_data_dir: the directory to annotation file
# - Save_data_dir: the directory to save the filtered data
# Output
# - kegg_level3_pathway_all.csv: all CO2 fixation related KEGG level 3 pathways
# - MAG_CO2_fixation_gene_num_summary.csv: summary of MAGs with CO2 fixation gene number
# """

# %% 
import os
import pandas as pd

# 输入输出路径
Anno_data_dir = "../../data/Fig2/genome_annotation/"
Save_data_dir = "../../data/Fig2/kegg_level3_pathway/"

# ===== 新增参数：是否提取所有代谢通路，还是只提取特定通路 =====
# 设置 extract_all_pathways = True 提取所有代谢通路，否则只提取关心的通路
extract_all_pathways = True  # 修改为True则提取所有Level3通路

# 你关心的KEGG level 3通路（仅在extract_all_pathways=False时生效）
pathways = [
    'Carbon fixation pathways in prokaryotes',
    'Carbon fixation in photosynthetic organisms'
]

kegg_data_path = 'kegg/kegg_level3_stat.xls'

# %%
# 用于存储所有MAG的结果
all_results = []

# 遍历所有以MAG开头的文件夹
for mag_folder in os.listdir(Anno_data_dir):
    mag_path = os.path.join(Anno_data_dir, mag_folder)
    if os.path.isdir(mag_path) and mag_folder.startswith("MAG"):
        kegg_file = os.path.join(mag_path, kegg_data_path)
        if os.path.exists(kegg_file):
            try:
                df = pd.read_csv(kegg_file, sep='\t', dtype=str)
                # 根据参数决定是否提取所有Level3通路
                if extract_all_pathways:
                    filtered = df.copy()
                else:
                    filtered = df[df['Level3'].isin(pathways)].copy()
                if not filtered.empty:
                    filtered['MAG'] = mag_folder  # 添加MAG信息
                    all_results.append(filtered)
            except Exception as e:
                print(f"Error reading {kegg_file}: {e}")

# %%
# 合并所有结果
if all_results:
    merged_df = pd.concat(all_results, ignore_index=True)
    # 检查保存目录是否存在，不存在则创建
    if not os.path.exists(Save_data_dir):
        os.makedirs(Save_data_dir)
    # 保存为csv
    output_file = os.path.join(Save_data_dir, "kegg_level3_pathway_all.csv")
    merged_df.to_csv(output_file, index=False)
    print(f"合并结果已保存到: {output_file}")
else:
    print("未找到任何符合条件的KEGG注释文件。")

# %%
# 统计所有MAG中CO2固定相关（目标pathway）的基因数目并输出表格
# 仅在extract_all_pathways=False时统计CO2固定相关基因数
if all_results and not extract_all_pathways:
    merged_df = pd.concat(all_results, ignore_index=True)
    # 假设每一行代表一个基因，统计每个MAG的基因数
    co2_gene_counts_df = merged_df.groupby('MAG').size().reset_index(name='CO2_fixation_gene_num')
    co2_gene_counts_output = os.path.join(Save_data_dir, "MAG_CO2_fixation_gene_num_summary.csv")
    co2_gene_counts_df.to_csv(co2_gene_counts_output, index=False)
    print(f"所有MAG的CO2固定相关基因数目统计已保存到: {co2_gene_counts_output}")
elif not all_results:
    print("未找到任何MAG的CO2固定相关基因注释，无法统计基因数目。")
elif extract_all_pathways:
    print("extract_all_pathways=True，未统计CO2固定相关基因数目。")

# %%
