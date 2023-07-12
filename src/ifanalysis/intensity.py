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



def cellcycle_intensity(df: pd.DataFrame, conditions: list, intensity_col: str, label, title: str, colors: list = colors, path=path):
    """Function to plot combined violin plots and barplot of a specific measurment in G1, S and G2/M phases

    Args:
        df : normalised dataframe from cell_cycle_normalisation function
        conditions : list of conditions to plot
        intensity_col : meaure from IF experiment to plot
        label : name of the meauremnts fior the y axis and title
        title : Defaults to '{intensity} intensity in G1, S and G2/M phases'.
        colors : Defaults to standard qualitative palette
    """
    cellcycle_phases = ['G1', 'S', 'G2/M']
    fig, ax = plt.subplots(nrows= 3, figsize=(len(conditions), 5))
    for i, phase in enumerate(cellcycle_phases):
        df_phase = df[df['cell_cycle'] == phase]
        max_value = df_phase[intensity_col].quantile(0.999)
        sns.violinplot(x="condition", y=intensity_col, data=df_phase, order=conditions, dodge=False, scale="width", palette=[colors[i]], inner=None, linewidth=0, ax=ax[i])
        for violin in ax[i].collections:
            violin.set_alpha(1)
        sns.boxplot(x="condition", y=intensity_col, data=df_phase, order=conditions, showfliers=False, color='white', saturation=0.5, width=0.2, boxprops={'zorder': 2},ax=ax[i])
        ax[i].set_ylabel(f"{label} {phase}")
        if i != 2:
            ax[i].set_xticklabels([])
        else:
            ax[i].set_xticklabels(conditions, rotation=30, ha='right')
        ax[i].set_xlabel('') 
        ax[i].set_ylim(0, max_value)
        fig.suptitle(title, size=16, weight='bold', x=0.05, horizontalalignment='left')
        fig._suptitle.set_weight('bold')
    save_fig(path, title)