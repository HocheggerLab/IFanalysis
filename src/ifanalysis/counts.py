import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
from ifanalysis.normalisation import *
from ifanalysis._helper_functions import save_fig
from typing import List, Tuple

# Matplotlib Style and Colors

STYLE = Path('~/matplotlib_style/Style_01.mplstyle')
plt.style.use(STYLE)
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]

# Get the path to the styles directory within the package
path = Path.cwd()

def rel_cellnumber(df_count: pd.DataFrame, conditions: List[str], cell_line: str, colors: List[str] = colors,
                   title: str = None, save: bool = False, path: Path =path, ax=None):
    """
    Function to plot relative cell number per condition and cell line
    :param df_count: dataframe provided by count_per_cond function,
    :param conditions: List of conditions to be plotted
    :param cell_line: cell line to be plotted
    :param colors: option, List of colors to be used for plotting. default, colors from prop_cycle
    :param title: option, name of the figure used as title in matplotlib and to save the figure default: 'Relative Cell Number'
    :param save: boolean, default True, if True saves the figure in the path provided
    :param path: option, default: Path.cwd(), path to save the figure
    :return: None, saves figure in path
    """
    ax = ax or plt.gca()
    plot = sns.barplot(
        data= df_count[df_count.cell_line == cell_line],
        x="condition",
        y="norm_count",
        errorbar="sd",
        order=conditions,
        palette=colors,
        ax=ax,
    )
    if title:
        plt.title(f"{title} {cell_line}")
    if save:
        save_fig(path, f"{title} {cell_line}")
    return plot

def abs_cellnumber(df_count: pd.DataFrame, conditions, cell_line: str, colors = colors,
                   title: str ='Absolute Cell Number', save: bool = True, path: Path = path) -> None:
    """
    Function to plot absulute cell number per condition and cell line
    :param df_count: dataframe provided by count_per_cond function,
    :param conditions: List of conditions to be plotted
    :param cell_line: cell line to be plotted
    :param colors: option, List of colors to be used for plotting. default, colors from prop_cycle
    :param title: option, name of the figure used as title in matplotlib and to save the figure default: 'Relative Cell Number'
    :param save: boolean, default True, if True saves the figure in the path provided
    :param path: option, default: Path.cwd(), path to save the figure
    :return: None, saves figure in path
    """
    fig, ax = plt.subplots(figsize=(len(conditions), 3))
    sns.barplot(
        data= df_count[df_count.cell_line == cell_line],
        x="condition",
        y="abs cell count",
        errorbar="sd",
        order=conditions,
        palette=colors,
        ax=ax,
    )
    ax.set_title(f"{title} {cell_line}")
    if save:
        save_fig(path, f"{title} {cell_line}")

def count_plots(df_count: pd.DataFrame, conditions: List[str]):
    row_num = len(df_count.cell_line.unique())
    fig, ax = plt.subplots(nrows=row_num, figsize=(6, 2*row_num))
    if row_num == 1:
        ax = [ax]  # Wrap single
    for cell_line, ax_single in zip(df_count.cell_line.unique(), ax):
        rel_cellnumber(df_count, conditions, cell_line, ax=ax_single)
        ax_single.set_title(cell_line)
        ax_single.set_xlabel('')
    plt.tight_layout()
