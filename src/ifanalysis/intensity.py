import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
from ifanalysis.normalisation import *
from ifanalysis._helper_functions import save_fig
from typing import List, Tuple

# Matplotlib Style and Colors
# Get the path to the styles directory within the package
STYLE = Path('~/matplotlib_style/Style_01.mplstyle')
plt.style.use(STYLE)
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]


path = Path.cwd()



def plot_int_scatter(ax, df, condition, col, cell_number): 
    phases = ["Sub-G1", "G1", "S", "G2/M", "Polyploid"]
    data = df[df['condition'] == condition]
    if cell_number and len(data) >= cell_number:
            data = data.sample(n=cell_number, random_state=42)
    sns.scatterplot(
        data=data, 
        x="integrated_int_DAPI_norm", 
        y=col, 
        hue='cell_cycle',
        hue_order=phases,
        s=5,
        alpha=0.8,
        ax=ax)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int(x))))
    ax.set_xticks([1, 2, 4, 8])
    ax.set_xlim([1, 16])
    ax.set_xlabel("norm. DNA content")
    ax.set_title(condition)
    ax.legend().remove()

def plot_int_violin(ax, i, df: pd.DataFrame, conditions: list, col: str, phase):     
    sns.violinplot(
        x="condition", 
        y=col, 
        data=df, 
        order=conditions, 
        dodge=False, 
        density_norm='width', 
        inner=None, 
        linewidth=0,
        color=colors[i+1],
        ax=ax)
    for violin in ax.collections:
        violin.set_alpha(1)
    sns.boxplot(x="condition", y=col, data=df, order=conditions, showfliers=False, color='white', saturation=0.5, width=0.2, boxprops={'zorder': 2},ax=ax)
    # Set the tick positions and labels before customizing the tick labels
    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels(conditions, rotation=30, ha='right')
    
    ax.set_ylabel('')
    ax.set_title(phase)   
    ax.set_xlabel('') 

def int_combplot(df, conditions, cell_line, col, label=None, title=None, cellnumber=None, save=True, path=path):
    
    col_number = len(conditions)
    if label not in col:
        df1 = df[(df.cell_line == cell_line) & (df.label == label)]
    else:
        df1 = df[(df.cell_line == cell_line)]
    condition_list = conditions * 2
    y_max = df1[col].quantile(0.99) * 1.5
    y_min = df1[col].quantile(0.01) * 0.8
    fig = plt.figure(figsize=(3 * len(conditions), 3))
    gs = GridSpec(1, col_number + 3, width_ratios=[1] * col_number + [0.5] * 3, wspace=0.1)
    for i, condition in enumerate(conditions):
        ax_int = fig.add_subplot(gs[0, i])
        ax_int.set_ylim([y_min+1, y_max+1])
        plot_int_scatter(ax_int, df1, condition, col, cellnumber)
        if i == 0:
            ax_int.set_ylabel(col)
            ax_int.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int(x))))
        else:    
            ax_int.set_ylabel('')
            ax_int.set_yticklabels([])
            ax_int.legend().remove()
    max_value = df1[col].quantile(0.995)
    min_value = 0
    cellcycle_phases = ['G1', 'S', 'G2/M']
    for i, phase in enumerate(cellcycle_phases):
        ax_violin = fig.add_subplot(gs[0, i+col_number])
        ax_violin.set_ylim([y_min, y_max])
        ax_violin.set_yticklabels([])
        df_phase = df1[df1['cell_cycle'] == phase]
        plot_int_violin(ax_violin, i, df_phase, conditions, col, phase)
    if not title:
        title = f"{col}, {cell_line}"
    fig.suptitle(title, size=14, weight='bold', x=0.01, horizontalalignment='left')
    fig._suptitle.set_weight('bold')
    plt.subplots_adjust(top=0.8)  # Increase the value to leave more space
    if save:
        save_fig(fig, path, title, tight_layout=False, fig_extension="png")