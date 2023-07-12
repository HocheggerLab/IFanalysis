from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import seaborn as sns

from ifanalysis._helper_functions import save_fig
from ifanalysis.normalisation import *

# Matplotlib Style and Colors

STYLE = Path('~/matplotlib_style/Style_01.mplstyle')
plt.style.use(STYLE)
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]

# Get the path to the styles directory within the package
path = Path.cwd()


def count_per_cond(df: pd.DataFrame, ctr_cond: str) -> pd.DataFrame:
    """
    Function to generate counts per condition and cell line data using groupby
    of the single cell dataframe from omero-screen. Data are grouped by
    cell line and condition. The function also normalises
    the data using normalise count function with the supplied ctr_cond as a reference.
    :param df: dataframe from omero-screen
    :param ctr_cond: control condition
    :return: dataframe with counts per condition and cell line
    """
    df_count = (
        df.groupby(["cell_line", "condition", "well", "plate_id"])["experiment"]
        .count()
        .reset_index()
        .rename(columns={"experiment": "abs cell count"})
    )
    df_count["norm_count"] = normalise_count(df_count, ctr_cond)
    return df_count


def normalise_count(df: pd.DataFrame, ctr_cond: str) -> pd.Series:
    """
    Function to normalise counts per condition and cell line data using groupby
    :param df: grouped by dataframe provided by norm_count function
    :param ctr_cond: control condition
    :return: pandas series with normalised counts
    """
    norm_count = pd.DataFrame()
    for cell_line in df["cell_line"].unique():
        norm_value = df.loc[
            (df["cell_line"] == cell_line) & (df["condition"] == ctr_cond),
            "abs cell count",
        ].mean()
        rel_cellcount = (
                df.loc[df["cell_line"] == cell_line, "abs cell count"] / norm_value
        )
        norm_count = pd.concat([norm_count, rel_cellcount])
    return norm_count


def cellnumber(df_count: pd.DataFrame, conditions: List[str], cell_line: str, type: str = "relative",
                   title: str = None, save: bool = False, path: Path = None, ax=None):
    """
    Function to plot relative cell number per condition and cell line

    :param df_count: dataframe provided by count_per_cond function,
    :param conditions: List of conditions to be plotted
    :param cell_line: cell line to be plotted
    :param title: option, default: None
    :param save: boolean, default False, if True saves the figure in the path provided
    :param path: option, default: None, path to save the figure
    :return: matplotlib axis object
    """
    if type == "relative":
        y_axis = "norm_count"
    elif type == "absolute":
        y_axis = "abs cell count"
    else:
        raise ValueError("type must be relative or absolute")
    ax = ax or plt.gca()
    plot = sns.barplot(
        data=df_count[df_count.cell_line == cell_line],
        x="condition",
        y=y_axis,
        errorbar="sd",
        order=conditions,
        ax=ax,
    )
    if title:
        plt.title(f"{title} {cell_line}")
    if save:
        save_fig(path, f"{title} {cell_line}")
    return plot



def count_plots(df_count: pd.DataFrame, conditions: List[str], title, save: bool = True, path: Path = path):
    """
    Function to plot relative and absolute cell number per condition and cell line
    combined in one figure
    :param df_count: pandas df from count per condition function
    :param conditions: list of conditions to be plotted
    :param save: boolean, default True, if True saves the figure in the path provided
    :param path: option, default: None, path to save the figure
    :return: None (plots and saves the figure)
    """
    row_num = len(df_count.cell_line.unique())
    fig, ax = plt.subplots(nrows=row_num, ncols=2, figsize=(6, 2 * row_num))
    ax = np.atleast_2d(ax)  # Reshape ax to always have 2 dimensions

    for i, cell_line in enumerate(df_count.cell_line.unique()):
        cellnumber(df_count, conditions, cell_line, type='relative', ax=ax[i, 0])
        ax[i, 0].set_title(cell_line)
        ax[i, 0].set_xlabel('')
        ax[i, 0].set_xticklabels(conditions, rotation=30, ha='right')
        cellnumber(df_count, conditions, cell_line, type='absolute', ax=ax[i, 1])
        ax[i, 1].set_xlabel('')
        ax[i, 1].set_xticklabels(conditions, rotation=30, ha='right')
    fig.suptitle(title, size=16, weight='bold', x=0.05, horizontalalignment='left')
    fig._suptitle.set_weight('bold')
    plt.tight_layout()
    if save:
        save_fig(path, title)

def count_cells(df: pd.DataFrame, conditions: List[str], ctr_cond: str = 'CTR', title: str ='Cell Counts"',
                save: bool = True, path: Path = path):
    """
    Function to count cells per condition and cell line and generate plots
    :param df: dataframe from omero_screen
    :param conditions: list of conditions to be plotted
    path: path to save the figure
    :return:
    """
    df_count = count_per_cond(df, ctr_cond)
    count_plots(df_count, conditions, title, save=save, path=path)
