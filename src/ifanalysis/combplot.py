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



def plot_histogram(ax, i, data):
    sns.histplot(data=data, x="integrated_int_DAPI_norm", ax=ax)
    ax.set_xlabel(None)
    ax.set_xscale("log", base=2)
    ax.set_xlim([1, 16])
    ax.xaxis.set_visible(False)   
    if i == 0:
        ax.set_ylabel("Frequency")
    else:
        ax.yaxis.set_visible(False)

def plot_scatter(ax, i, data, conditions, H3:bool):    
    if H3:
        phases = ["Sub-G1", "G1", "S", "G2" "M", "Polyploid"]
    else:
        phases = ["Sub-G1", "G1", "S", "G2/M", "Polyploid"]
    sns.scatterplot(
        data=data, 
        x="integrated_int_DAPI_norm", 
        y='intensity_mean_EdU_nucleus_norm', 
        hue='cell_cycle',
        hue_order=phases, 
        s=5,
        alpha=0.8,
        ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log', base=2)
    ax.grid(False)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int(x))))
    ax.set_xticks([1, 2, 4, 8])
    ax.set_xlim([1, 16])
    ax.set_xlabel("norm. DNA content (log2)")
    if i == len(conditions):
        ax.set_ylabel("norm. EdU intesity (log)")
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int(x))))
        ax.legend(loc='upper right')
    else:
        ax.yaxis.set_visible(False)
        ax.legend().remove()
    ax.axvline(x=3, color='black', linestyle='--')
    ax.axhline(y=3, color='black', linestyle='--')
    sns.kdeplot(
        data=data, 
        x="integrated_int_DAPI_norm", 
        y='intensity_mean_EdU_nucleus_norm', 
        fill=True,
        alpha=0.3,
        cmap='rocket_r', 
        ax=ax
)
   

def combplot(df, conditions, cell_line, title_str, cell_number=3000, H3=False, save=True, path=path):
    
    col_number = len(conditions)
    df1 = df[df.cell_line == cell_line]
    condition_list = conditions * 2

    fig = plt.figure(figsize=(3 * len(conditions), 5))
    gs = GridSpec(2, col_number, height_ratios=[1, 3])
    ax_list = [(i, j) for i in range(2) for j in range(col_number)]
    y_max = df['intensity_mean_EdU_nucleus_norm'].quantile(0.99) * 1.5
    y_min = df['intensity_mean_EdU_nucleus_norm'].quantile(0.01) * 0.8

    for i, pos in enumerate(ax_list):
        data = df1[df1.condition == condition_list[i]]
        if len(data) >= cell_number:
            data_red = data.sample(n=cell_number, random_state=42)
        else:
            data_red = data
        ax = fig.add_subplot(gs[pos[0], pos[1]])

        if i < len(conditions):
            plot_histogram(ax, i, data_red)
            ax.set_title(f"{condition_list[i]}, {len(data_red)} cells", size=12, weight='bold')
        else:
            plot_scatter(ax, i, data_red, conditions, H3)
            ax.set_ylim([y_min, y_max])

        ax.grid(visible=False)
    title = f"{title_str} {cell_line}"
    fig.suptitle(title, size=16, weight='bold', x=0.05, horizontalalignment='left')
    fig._suptitle.set_weight('bold')
    if save:
        save_fig(path, f"CombPlot_{title_str}_{cell_line}", fig_extension="png")

