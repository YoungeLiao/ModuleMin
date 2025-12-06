# %%
import pandas as pd
import os
import matplotlib.pyplot as plt

def load_data(file_path, sheet_name=0):
    """Load data from an Excel file."""
    try:
        data = pd.read_excel(file_path, sheet_name=sheet_name)
        return data
    except ValueError as e:
        print(f"Error loading {file_path}: {e}")
        return None

def add_sequence_id_column(data):
    """Add a new column 'Sequence_ID' by removing the last two characters of the first column."""
    data['Sequence_ID'] = data.iloc[:, 0].str[:-2]
    return data

def merge_data(data1, data2, on_column):
    """Merge two DataFrames based on a specified column."""
    return pd.merge(data1, data2, on=on_column)

def save_data(data, output_file_path):
    """Save the DataFrame to an Excel file."""
    data.to_excel(output_file_path, index=False)
    print(f"Merged data has been saved to {output_file_path}")

def plot_pie_chart_top_10(data, genus_column, total_column, output_file_path):
    """绘制美化后的前10属饼图，其余合并为Other，符合Nature风格。去除属名前的g__前缀。"""
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    # 设置全局字体
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['font.size'] = 16
    mpl.rcParams['axes.linewidth'] = 1.2

    # 去除g__前缀，仅保留属名，并处理特殊属名
    data = data.copy()
    data[genus_column] = (
        data[genus_column]
        .astype(str)
        .str.replace(r'^g__', '', regex=True)
        .str.replace('Candidatus_Scalindua', 'Candidatus Scalindua')
    )

    genus_data = data.groupby(genus_column)[total_column].sum()
    top_10_genus = genus_data.nlargest(10)
    other_genus = genus_data[~genus_data.index.isin(top_10_genus.index)].sum()
    top_10_genus['Other'] = other_genus

    # 莫兰迪色系基础上提升亮度和饱和度，保证分组易区分
    colors = [
        "#7FC8C0",  # 明亮蓝绿
        "#F2B5A8",  # 明亮粉橙
        "#A8E6A3",  # 明亮嫩绿
        "#FFD6A6",  # 明亮浅橙
        "#A6D8FF",  # 明亮天蓝
        "#D1A6FF",  # 明亮紫色
        "#E6A6C9",  # 明亮粉紫
        "#A6FFE2",  # 明亮青色
        "#E6FFA6",  # 明亮浅绿
        "#FFB6B6",  # 明亮粉红
        "#FFF2D7"   # 明亮米色
    ]

    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
    wedges, texts = ax.pie(
        top_10_genus,
        labels=None,  # 不直接在饼图上显示标签，防止重叠
        autopct=None,
        startangle=90,
        colors=colors,
        pctdistance=0.75,
        labeldistance=1.05,
        wedgeprops=dict(linewidth=1.2, edgecolor='white')
    )
    # 计算每个扇区的角度和中心点坐标
    import numpy as np

    total = sum(top_10_genus)
    theta1 = 90  # 起始角度
    label_positions = []
    for v in top_10_genus:
        theta2 = theta1 + v / total * 360
        mid_angle = (theta1 + theta2) / 2
        label_positions.append(mid_angle)
        theta1 = theta2

    # 只允许一个Other标签，避免重复
    other_labeled = False
    for i, (label, value, angle) in enumerate(zip(top_10_genus.index, top_10_genus.values, label_positions)):
        percent = value / total * 100
        if percent < 2:
            continue  # 跳过小于2%的标签
        # 斜体属名，百分比不斜体
        genus_label = r"$\it{" + label.replace(' ', r'\ ') + "}$" if label != 'Other' else label
        if i < 4:
            # 直接贴在扇区外
            x = 0.85 * np.cos(np.deg2rad(angle))
            y = 0.85 * np.sin(np.deg2rad(angle))
            ha = 'left' if x > 0 else 'right'
            va = 'center'
            ax.text(
                x, y, f"{genus_label}\n{percent:.1f}%", ha=ha, va=va, fontsize=14,
                bbox=None  # 无背景
            )
        elif i == 4 and label != 'Other':
            # top5用短引线，线更靠近扇区（不对Other用引线）
            x = 0.95 * np.cos(np.deg2rad(angle))
            y = 0.95 * np.sin(np.deg2rad(angle))
            lx = 1.10 * np.cos(np.deg2rad(angle))
            ly = 1.10 * np.sin(np.deg2rad(angle))
            ha = 'left' if lx > 0 else 'right'
            va = 'center'
            ax.plot([x, lx], [y, ly], color='gray', lw=1)
            ax.text(
                lx, ly, f"{genus_label}\n{percent:.1f}%", ha=ha, va=va, fontsize=14,
                bbox=None  # 无背景
            )
        elif label == 'Other' and not other_labeled:
            # 只标注一次Other
            x = 0.85 * np.cos(np.deg2rad(angle))
            y = 0.85 * np.sin(np.deg2rad(angle))
            ha = 'left' if x > 0 else 'right'
            va = 'center'
            ax.text(
                x, y, f"{genus_label}\n{percent:.1f}%", ha=ha, va=va, fontsize=14,
                bbox=None  # 无背景
            )
            other_labeled = True
        # top6~10和重复的Other不标注

    # 去除边框和多余元素
    ax.set(aspect="equal")
    plt.tight_layout()
    plt.savefig(output_file_path, bbox_inches='tight', transparent=True)
    plt.show()


def plot_bar_top10_genus_by_sample(
    merged_data, 
    genus_column='Genus', 
    sample_columns=None, 
    metadata_path='../data/metadata.xlsx', 
    output_file_path='../results/barplot_top10_genus.pdf'
):
    """
    绘制前10属在不同组（如Dark和IR）中的丰度柱状图，横坐标为属，分组为颜色，带误差棒和数据点，风格与饼图一致。
    支持纵坐标丰度取对数，便于比较低丰度属。
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np

    # 设置全局字体和风格，保持与饼图一致
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['font.size'] = 16
    mpl.rcParams['axes.linewidth'] = 1.2

    if sample_columns is None:
        sample_columns = ['IR1', 'IR2', 'IR3', 'Dark1', 'Dark2', 'Dark3']

    # 获得前10属
    genus_sample = merged_data.groupby(genus_column)[sample_columns].sum()
    top10 = genus_sample.sum(axis=1).nlargest(10).index
    genus_sample_top10 = genus_sample.loc[top10]

    # 转为长表
    genus_sample_long = genus_sample_top10.reset_index().melt(
        id_vars=genus_column, value_vars=sample_columns,
        var_name='Sample', value_name='Abundance'
    )

    # 加载metadata，获取分组信息
    metadata = pd.read_excel(metadata_path)
    if 'Sample' not in metadata.columns:
        metadata.columns = [str(c) for c in metadata.columns]
        if 'sample' in metadata.columns:
            metadata.rename(columns={'sample': 'Sample'}, inplace=True)
    genus_sample_long = pd.merge(genus_sample_long, metadata, on='Sample', how='left')

    # 分组列
    group_col = 'Group' if 'Group' in metadata.columns else 'group'
    group_order = genus_sample_long[group_col].drop_duplicates().tolist()

    # 莫兰迪色系，和饼图一致
    colors = [
        "#d67875", "#818181","#7FC8C0","#A8E6A3", "#FFD6A6", 
        "#D1A6FF", "#E6A6C9", "#A6FFE2", "#E6FFA6", "#FFB6B6"
    ]
    palette = dict(zip(group_order, colors[:len(group_order)]))

    # 计算均值和标准差（对数前先加1防止log(0)）
    genus_sample_long['Abundance_log'] = np.log10(genus_sample_long['Abundance'] + 1)
    summary = genus_sample_long.groupby([genus_column, group_col])['Abundance_log'].agg(['mean', 'std']).reset_index()

    # 画图
    fig, ax = plt.subplots(figsize=(12, 6), dpi=300)
    width = 0.35  # 每组的宽度
    n_group = len(group_order)
    x = np.arange(len(top10))  # 属的横坐标

    for i, group in enumerate(group_order):
        group_data = summary[summary[group_col] == group]
        # 保证顺序
        group_data = group_data.set_index(genus_column).reindex(top10).reset_index()
        bar_x = x + (i - (n_group-1)/2) * width
        ax.bar(
            bar_x, group_data['mean'], width=width, 
            yerr=group_data['std'], label=group, color=palette[group], capsize=3, edgecolor='white', linewidth=1.2
        )
        # 叠加数据点
        for j, genus in enumerate(top10):
            points = genus_sample_long[(genus_sample_long[genus_column] == genus) & (genus_sample_long[group_col] == group)]['Abundance_log']
            # 抖动
            jitter = np.random.uniform(-width/3, width/3, size=len(points))
            ax.scatter(
                np.full(len(points), bar_x[j]) + jitter, points, 
                color=palette[group], edgecolor='k', s=30, alpha=0.7, zorder=10
            )

    ax.set_xticks(x)
    ax.set_xticklabels(top10, rotation=30, ha='right')
    ax.set_xlabel('')
    ax.set_ylabel('log Abundance')
    ax.set_title('IR Abundance in Top 10 Genus')
    # 修改legend位置到图内右上角
    ax.legend(title='Group', loc='upper right', frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_file_path, bbox_inches='tight', transparent=True)
    plt.show()

if __name__ == "__main__":
    main()

# # === DEBUGS ===
# # %%
# import pandas as pd

# def load_data(file_path):
#     """Load data from an Excel file."""
#     return pd.read_excel(file_path)

# annotated_ir_file = "../data/Annotated_IROptogenetics.xlsx"
# annotated_ir_data = load_data(annotated_ir_file)
# # %%

# %%
