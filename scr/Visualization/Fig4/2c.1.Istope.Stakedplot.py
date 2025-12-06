# %% 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D

def plot_isotopic_stacking(data_path, save_path, figsize=(3, 5)):
    """
    Creates Nature-style stacked bar plot of isotopic patterns

    Args:
        data_path (str): Path to Excel file with isotope percentages
        save_path (str): Output path for the figure
        figsize (tuple): Figure dimensions in inches (width, height)
        dpi (int): Resolution for output figure
    """
    # Nature-style configuration (font size per Nature guidelines)
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 12,
        'axes.labelsize': 12,
        'axes.titlesize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'figure.dpi': 300,
        'figure.figsize': (5.4/2.54, 4.5/2.54),  # 1.5倍单栏尺寸（宽13.5cm, 高11.25cm）
        'axes.linewidth': 1,
        'lines.linewidth': 1.2,
        'grid.linewidth': 0.4,
        'grid.alpha': 0.2,
        'savefig.dpi': 600,
        'savefig.format': 'pdf',
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.08
    })

    # Load and preprocess data
    df = pd.read_excel(data_path)

    # Prepare data matrix for plotting
    categories = ['Dark_AcCoA', 'IR_AcCoA', 'Dark_MalCoA', 'IR_MalCoA']
    data_matrix = np.array([
        df.loc[df['Sample'].str.contains('AcCoA') & df['Sample'].str.contains('Dark'), 'Unlabeled':'m3'].values.flatten(),
        df.loc[df['Sample'].str.contains('AcCoA') & df['Sample'].str.contains('IR'), 'Unlabeled':'m3'].values.flatten(),
        df.loc[df['Sample'].str.contains('MalCoA') & df['Sample'].str.contains('Dark'), 'Unlabeled':'m3'].values.flatten(),
        df.loc[df['Sample'].str.contains('MalCoA') & df['Sample'].str.contains('IR'), 'Unlabeled':'m3'].values.flatten()
    ])

    # Plotting
    fig, ax = plt.subplots(figsize=figsize)

    # Nature-style colors
    n_isotopes = data_matrix.shape[1]
    colors = ['#F5EDED', '#818FB4', '#435585', '#363062']
    hatches = ['', '', '', '//']
    if n_isotopes > len(colors):
        colors = (colors * ((n_isotopes // len(colors)) + 1))[:n_isotopes]
    if n_isotopes > len(hatches):
        hatches = (hatches * ((n_isotopes // len(hatches)) + 1))[:n_isotopes]

    # Stacked bar plot
    bottom = np.zeros(len(categories))
    for i in range(data_matrix.shape[1]):
        bars = ax.bar(categories, data_matrix[:, i], bottom=bottom,
                      color=colors[i], edgecolor='white', hatch=hatches[i])
        bottom += data_matrix[:, i]

    # Axis and title formatting
    ax.set_ylim(0, 100.5)
    ax.set_ylabel('Isotopic Abundance (%)', labelpad=8)
    ax.set_title('', pad=10)
    ax.yaxis.grid(True, linestyle='--', alpha=0.4)

    # x轴文字旋转避免重合
    ax.set_xticklabels(categories, rotation=30, ha='right')

    # # Condition labels（放到下方）
    # ax.text(-0.5, -13, 'DARK', fontsize=8, ha='center', fontweight='bold', color='#1f77b4')
    # ax.text(0.5, -13, 'IR', fontsize=8, ha='center', fontweight='bold', color='#d62728')
    # ax.text(1.5, -13, 'DARK', fontsize=8, ha='center', fontweight='bold', color='#1f77b4')
    # ax.text(2.5, -13, 'IR', fontsize=8, ha='center', fontweight='bold', color='#d62728')

    # Add significance indicators
    ax.plot([0, 1], [102, 102], color='k', linewidth=1)
    ax.plot([2, 3], [102, 102], color='k', linewidth=1)

    # Custom legend（无边框）
    legend_elements = [
        Line2D([0], [0], color='w', marker='s', markersize=8, markerfacecolor=colors[0], label='Unlabeled'),
        Line2D([0], [0], color='w', marker='s', markersize=8, markerfacecolor=colors[1], label='m1'),
        Line2D([0], [0], color='w', marker='s', markersize=8, markerfacecolor=colors[2], label='m2'),
        Line2D([0], [0], color='w', marker='s', markersize=8, markerfacecolor=colors[3], label='m3')
    ]
    ax.legend(handles=legend_elements, title='',
              frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=4)

    # 调整x轴下边距，避免下方标签被裁剪
    plt.subplots_adjust(bottom=0.22)

    # Save figure
    plt.savefig(save_path)
    print(f'Isotopic patterns saved to {save_path}')

if __name__ == '__main__':
    # Configurable parameters
    DATA_PATH = '../../../data/Fig4/e_Isotope.percentage.xlsx'
    SAVE_PATH = '../results/isotopic_stacking.pdf'

    plot_isotopic_stacking(DATA_PATH, SAVE_PATH)

# %%
