import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import ticker
import seaborn as sns
from ifanalysis._helper_functions import save_fig
from pathlib import Path
from typing import Optional
path = Path.cwd() 


def plot_histogram(ax, i, data):
    sns.histplot(data=data, x="integrated_int_DAPI_norm", ax=ax)
    ax.set_xlabel(None)
    ax.set_xscale("log", base=2)
    ax.set_xlim([1, 16])
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(int(x))))
    ax.set_xticks([1, 2, 4, 8])
    ax.set_xlabel("norm. DNA content")
    if i == 0:
        ax.set_ylabel("Frequency")
    else:
        ax.yaxis.set_visible(False)

def combplot_hist(df, conditions, title_str, cell_number=None, save=True, path=path):
    
    col_number = len(conditions)

    fig = plt.figure(figsize=(2 * len(conditions), 3))
    gs = GridSpec(1, col_number)  # Only one row of subplots

    for i, condition in enumerate(conditions):
        data = df[df.condition == condition]
        if cell_number and len(data) >= cell_number:
            data_red = data.sample(n=cell_number, random_state=42)
        else:
            data_red = data 
        ax = fig.add_subplot(gs[0, i])  # gs indexing changed to single row

        plot_histogram(ax, i, data_red)
        ax.set_title(f"{condition} \n{len(data_red)} cells", size=12, weight='bold')
        ax.grid(visible=False)
    
    fig.suptitle(title_str, size=14, weight='bold', x=0.01, horizontalalignment='left')
    fig._suptitle.set_weight('bold')
    plt.subplots_adjust(top=0.85)
    plt.tight_layout()
    if save:
        save_fig(fig, path, title_str, fig_extension="png")

def plot_int_violin(ax, df: pd.DataFrame, conditions: list, col: str, hue: Optional[str]):
    # Determine the number of hue levels

    if hue:
        sns.violinplot(
            x="condition", 
            y=col, 
            data=df, 
            order=conditions,
            dodge=True, 
            density_norm='width', 
            inner=None, 
            linewidth=0,
            hue=hue,
            ax=ax)
        
         # Calculate the width of each violin plot portion
        hue_levels = df[hue].unique()
        hue_count = len(hue_levels)
        width_per_hue = 0.8 / hue_count  # 0.8 is the default width of violin plots in seaborn

        for i, condition in enumerate(conditions):
            for j, level in enumerate(hue_levels):
                subset = df[(df["condition"] == condition) & (df[hue] == level)]
                # Calculate the center position for the boxplot
                box_center = i - 0.4 + width_per_hue/2 + j*width_per_hue
                sns.boxplot(
                    y=subset[col], 
                    showfliers=False, 
                    color='white', 
                    saturation=0.5, 
                    width=0.1, 
                    boxprops={'zorder': 2},
                    positions=[box_center], 
                    ax=ax)
    else:
        sns.violinplot(
            x="condition", 
            y=col, 
            data=df, 
            order=conditions,
            dodge=False, 
            density_norm='width', 
            inner=None, 
            linewidth=0,
            ax=ax)
        sns.boxplot(
            x="condition", 
            y=col, 
            data=df, 
            order=conditions, 
            showfliers=False, 
            color='white', 
            saturation=0.5, 
            width=0.2, 
            boxprops={'zorder': 2},
            ax=ax)

    for violin in ax.collections:
        violin.set_alpha(1)

    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels(conditions, rotation=30, ha='right')
    ax.set_ylabel('')
    ax.set_xlabel('')


def intensity_plot(df, conditions, columns, hue=None, title=None, save=True, path=path):
    fig, axarr = plt.subplots(nrows=len(columns), figsize=(len(conditions), len(columns)*3))
    
    # Ensure axarr is always a list
    if len(columns) == 1:
        axarr = [axarr]
    
    for i, col in enumerate(columns):
        y_max = df[col].quantile(0.99) * 1.5
        y_min = df[col].quantile(0.01) * 0.8
        plot_int_violin(axarr[i], df, conditions, col, hue)
        axarr[i].set_ylim([y_min, y_max])
        axarr[i].set_title(col)   

    fig.suptitle(title, size=14, weight='bold', x=0.01, horizontalalignment='left')
    plt.subplots_adjust(top=0.8)
    if save:
        save_fig(fig, path, title, tight_layout=True)
    