# %%
import pandas as pd
import seaborn as sns
import numpy as np

# 读取数据
asv = pd.read_excel('../../data/rawdata.shared/asv_taxon_all.xlsx')
mag = pd.read_excel('../../data/rawdata.shared/all.MAG.summary.xlsx')
mag_abundance = pd.read_excel('../../data/rawdata.shared/MAGs_rpkm_abundance.xlsx')
ir = pd.read_excel('../../results_data/Fig2/Database/Annotated_IROptogenetics.xlsx')

# 1. 统计IR基因元件中Genus的基因数（Gene Counts，前top5）
ir['Genus'] = ir['Genus'].str.replace('g__', '', regex=False)
ir_genus_counts = ir.groupby('Genus').size().sort_values(ascending=False)
top_genus = ir_genus_counts.head(5).index.tolist()

# 2. 统计16S中这些Genus的基因数（即ASV条目数）
asv['Genus'] = asv['Genus'].str.replace('g__', '', regex=False)
# Anaerolinea特殊处理
asv.loc[asv['Genus'] == 'unclassified_c__Anaerolineae', 'Genus'] = 'Anaerolinea'
asv_genus_counts = asv[asv['Genus'].isin(top_genus)].groupby('Genus').size()

# 3. 处理MAG的Genus/Class注释
# 合并mag和mag_abundance，确保MAG丰度信息带有taxonomy注释
mag_abundance = mag_abundance.rename(columns={'MAG': 'Bin'})  # 确保MAG名列一致
mag = mag.rename(columns={'MAG': 'Bin'})  # 统一MAG名列
mag_merged = pd.merge(mag_abundance, mag[['Bin', 'Genus', 'Class']], on='Bin', how='left')

mag_merged['Genus'] = mag_merged['Genus'].str.replace('g__', '', regex=False)
mag_merged['Class'] = mag_merged['Class'].str.replace('c__', '', regex=False)
# 特殊处理
mag_merged['Genus'] = mag_merged['Genus'].replace({'SCAELEC01': 'Candidatus_Scalindua'})
mag_merged.loc[mag_merged['Class'] == 'Anaerolineae', 'Genus'] = 'Anaerolinea'
# mag_merged['Total'] = mag_merged.iloc[:, 1:mag_merged.columns.get_loc('Domain')].sum(axis=1)
mag_abundance = mag_merged  # 后续统一用mag_abundance变量
# %%
# 读取metadata
meta = pd.read_excel('../../data/rawdata.shared/metadata.xlsx')

# 获取分组信息
group_col = 'Group'  # 假设metadata中分组列名为'Group'
sample_groups = meta.set_index('Sample')[group_col].to_dict()

# 处理IR数据：按Genus和分组统计丰度
ir_grouped = ir.copy()
ir_grouped = ir_grouped[ir_grouped['Genus'].isin(top_genus)]
# 提取样本列（即metadata中有分组信息的列）
sample_cols = [col for col in ir_grouped.columns if col in sample_groups]
# 将宽表转为长表
ir_melt = ir_grouped.melt(id_vars=['Genus'], value_vars=sample_cols, var_name='Sample', value_name='Total')
ir_melt['Group'] = ir_melt['Sample'].map(sample_groups)
ir_group_sum = ir_melt.groupby(['Genus', 'Group'])['Total'].sum().unstack(fill_value=0)

# %%
# 处理ASV数据：按Genus和分组统计丰度
asv_grouped = asv.copy()
asv_grouped = asv_grouped[asv_grouped['Genus'].isin(top_genus)]
# 提取样本列
sample_cols = [col for col in asv_grouped.columns if col in sample_groups]
asv_melt = asv_grouped.melt(id_vars=['Genus'], value_vars=sample_cols, var_name='Sample', value_name='Abundance')
asv_melt['Group'] = asv_melt['Sample'].map(sample_groups)
asv_group_sum = asv_melt.groupby(['Genus', 'Group'])['Abundance'].sum().unstack(fill_value=0)

# 处理MAG数据：按Genus和分组统计丰度
mag_grouped = mag_abundance.copy()
mag_grouped = mag_grouped[mag_grouped['Genus'].isin(top_genus)]
mag_sample_cols = [col for col in mag_grouped.columns if col in sample_groups]
mag_melt = mag_grouped.melt(id_vars=['Genus'], value_vars=mag_sample_cols, var_name='Sample', value_name='Abundance')
mag_melt['Group'] = mag_melt['Sample'].map(sample_groups)
mag_group_sum = mag_melt.groupby(['Genus', 'Group'])['Abundance'].sum().unstack(fill_value=0)

# 检查MAG部分缺失的原因
# 问题出在 mag_sample_cols 为空，导致 mag_melt 为空
# 打印 mag_grouped.columns 和 sample_groups.keys() 检查列名是否一致
print("mag_grouped.columns:", mag_grouped.columns.tolist())
print("sample_groups.keys():", list(sample_groups.keys()))

# 修正：有可能 mag_grouped 的样本列名与 metadata 不一致（如大小写、前后空格等）
# 统一列名格式
mag_grouped.columns = mag_grouped.columns.str.strip()
sample_groups_fixed = {k.strip(): v for k, v in sample_groups.items()}

# 重新提取样本列
mag_sample_cols = [col for col in mag_grouped.columns if col in sample_groups_fixed]
print("After strip, mag_sample_cols:", mag_sample_cols)

# 重新 melt
mag_melt = mag_grouped.melt(id_vars=['Genus'], value_vars=mag_sample_cols, var_name='Sample', value_name='Abundance')
mag_melt['Group'] = mag_melt['Sample'].map(sample_groups_fixed)
mag_group_sum = mag_melt.groupby(['Genus', 'Group'])['Abundance'].sum().unstack(fill_value=0)

# 合并结果，包含每个样本的丰度
result_list = []
methods = [
    ('Metagenomics', ir_melt, 'Total'),  # IR数据为宏基因组
    ('16S', asv_melt, 'Abundance'),
    ('MAG', mag_melt, 'Abundance')
]

for method, df, value_col in methods:
    for genus in top_genus:
        sub = df[df['Genus'] == genus]
        for sample in sub['Sample'].unique():
            group = sub[sub['Sample'] == sample]['Group'].iloc[0]
            val = sub[sub['Sample'] == sample][value_col].sum()
            result_list.append({
                'Genus': genus,
                'Method': method,
                'Group': group,
                'Sample': sample,
                'Abundance': val
            })

result_df = pd.DataFrame(result_list)

# 检查MAG部分是否缺失，若缺失则排查mag_melt和mag_grouped
if result_df[result_df['Method'] == 'MAG'].empty:
    print("Warning: MAG abundance results are still missing. Please check mag_melt and mag_grouped for issues.")
    print("mag_grouped shape:", mag_grouped.shape)
    print("mag_melt shape:", mag_melt.shape)
    print("mag_melt head:\n", mag_melt.head())
    print("top_genus:", top_genus)
    print("mag_grouped['Genus'] unique:", mag_grouped['Genus'].unique())
    print("mag_melt['Genus'] unique:", mag_melt['Genus'].unique())

# 保存详细结果
result_df.to_excel('All_IR_Genus_Abundance_BySample.xlsx', index=False)
print(result_df)

# %%
