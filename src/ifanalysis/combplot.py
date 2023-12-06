import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.gridspec import GridSpec
from matplotlib.legend_handler import HandlerTuple
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
        phases = ["Sub-G1", "G1", "S", "G2", "M", "Polyploid"]
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
    ax.set_xlabel("norm. DNA content")
    if i == len(conditions):
        ax.set_ylabel("norm. EdU intesity")
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int(x))))
    
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

def cellcycle_barplot(ax, df, conditions, H3): 
    df_mean, df_std = prop_pivot(df, conditions, H3)
    df_mean.plot(kind="bar", stacked=True, yerr=df_std, width=0.75, ax=ax)
    ax.set_ylim(0, 110)
    ax.set_xticklabels(conditions, rotation=30, ha='right')
    ax.set_xlabel('')  # Remove the x-axis label)
    if H3:
        legend = ax.legend(['Sub-G1', 'G1', 'S', 'G2', 'M', 'Polyploid'], title='CellCyclePhase')
        ax.set_ylabel('% of population')
    else:
        legend = ax.legend(['Sub-G1', 'G1', 'S', 'G2/M', 'Polyploid'], title='CellCyclePhase')
    # Get current handles and labels
    handles, labels = ax.get_legend_handles_labels()
    handles, labels = handles[::-1], labels[::-1]
    # Clear the current legend
    legend.remove()
    # Create a new legend with the reversed handles and labels
    legend = ax.legend(handles, labels, title='CellCyclePhase', loc='lower left')
    frame = legend.get_frame()
    frame.set_alpha(0.5)
    ax.set_ylabel('% of population')
    ax.grid(False)

def combplot(df, conditions, cell_line, title_str, cell_number=None, H3=False, save=True, path=path):
    
    col_number = len(conditions)
    df1 = df[df.cell_line == cell_line]
    condition_list = conditions * 2

    fig = plt.figure(figsize=(3 * len(conditions), 5))
    gs = GridSpec(2, col_number+1, height_ratios=[1, 3])
    ax_list = [(i, j) for i in range(2) for j in range(col_number)]
    y_max = df['intensity_mean_EdU_nucleus_norm'].quantile(0.99) * 1.5
    y_min = df['intensity_mean_EdU_nucleus_norm'].quantile(0.01) * 0.8

    for i, pos in enumerate(ax_list):
        data = df1[df1.condition == condition_list[i]]
        if cell_number and len(data) >= cell_number:
            data_red = data.sample(n=cell_number, random_state=42)
        else:
            data_red = data 
        ax = fig.add_subplot(gs[pos[0], pos[1]])

        if i < len(conditions):
            plot_histogram(ax, i, data_red)
            ax.set_title(f"{condition_list[i]} \n{len(data_red)} cells", size=12, weight='bold')
        else:
            plot_scatter(ax, i, data_red, conditions, H3)
            ax.set_ylim([y_min, y_max])

        ax.grid(visible=False)
    # Add the subplot spanning both rows
    ax_last = fig.add_subplot(gs[:, -1])
    ax_last.grid(visible=False)
    # Add the barplot to the subplot
    cellcycle_barplot(ax_last, df1, conditions, H3)
    # ax_last.set_title("Cell Cycle Dist.", size=12, weight='bold')

    title = f"{title_str} {cell_line}"
    fig.suptitle(title, size=14, weight='bold', x=0.01, horizontalalignment='left')
    fig._suptitle.set_weight('bold')
    # Adjust the spacing between the title and the plot
    plt.subplots_adjust(top=0.85)  # Increase the value to leave more space
    plt.tight_layout()
    if save:
        save_fig(fig, path, f"{cell_line} CombPlot {title_str}", fig_extension="png")

