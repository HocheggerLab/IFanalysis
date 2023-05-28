import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
from ifanalysis.normalisation import *
from typing import List, Tuple

# Matplotlib Style and Colors
plt.style.use("/Users/hh65/matplotlib_style/Style_01.mplstyle")
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]

path = Path.cwd()
def save_fig(
    path: Path, fig_id: str, tight_layout : bool = True, fig_extension: str = "pdf",
        resolution: int = 300) -> None:
    """
    coherent saving of matplotlib figures as pdfs (default)
    :param path: path for saving
    :param fig_id: name of saved figure
    :param tight_layout: option, default True
    :param fig_extension: option, default pdf
    :param resolution: option, default 300dpi
    :return: None, saves Figure in poth
    """

    dest = path / f"{fig_id}.{fig_extension}"
    print("Saving figure", fig_id)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(dest, format=fig_extension, dpi=resolution)

def rel_cellnumber(df_count: pd.DataFrame, conditions: List[str], cell_line: str, colors: List[str] = colors,
                   title: str ='Relative Cell Number', save: bool = True, path: Path =path) -> None:
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
    fig, ax = plt.subplots(figsize=(len(conditions), 3))
    sns.barplot(
        data= df_count[df_count.cell_line == cell_line],
        x="condition",
        y="norm_count",
        errorbar="sd",
        order=conditions,
        palette=colors,
        ax=ax,
    )
    ax.set_title(f"{title} {cell_line}")
    if save:
        save_fig(path, f"{title} {cell_line}")

def cellcycleplot_comb(df_cc: pd.DataFrame, conditions: List[str], cell_line: str, bins=1000,
                        title: str = 'Combined Cell Cycle Plot', colors: List[str]= colors,
                        col: str ='intensity_mean_EdU_nucleus_norm', save: bool =True, path: Path = Path.cwd()) -> None:
    """
    Function to plot combined cell cycle plot (Histgram and Scatterplot) for a given condition and cell line
    :param df_cc: df from cell_cycle_analysis function
    :param conditions: list of conditions to be plotted
    :param cell_line: cell line to be analysed
    :param bins: option, default 1000, number of bins for the histogram
    :param title: default 'Combined Cell Cycle Plot'
    :param colors: option, default colors from prop_cycle
    :param col: column to plot at the y-axis, default 'intensity_mean_EdU_nucleus_norm'
    :param save: boolean, default True, if True saves the figure in the path provided
    :param path: option, default: Path.cwd(), path to save the figure
    :return: None, saves figure in path
    """
    phases = ["G1", "S", "G2/M", "Polyploid", "Sub-G1"]
    col_number = len(conditions)
    df_cc = df_cc[df_cc.cell_line == cell_line]
    condition_list = conditions * 2
    fig = plt.figure(figsize=(3*len(conditions), 5))
    # define the grid layout with different height ratios
    gs = GridSpec(2, col_number, height_ratios=[1, 3])
    ax_list = [(i,j) for i in range(2) for j in range(col_number)]
    y_max = df_cc[col].quantile(0.99) * 1.5
    y_min = df_cc[col].quantile(0.01) * 0.8
    for i, pos in enumerate(ax_list):
        data = df_cc[df_cc.condition == condition_list[i]]
        ax = fig.add_subplot(gs[pos[0], pos[1]])
        if i < len(condition_list)/2:
            data['integrated_int_DAPI_norm'].plot.hist(bins=bins, color=colors[-1], ax=ax)
            ax.set_title(condition_list[i])
            ax.set_xlabel(None)
            if i == 0:
                ax.set_ylabel("Frequency")
            else:
                ax.set_ylabel('')
        else:
            ax.scatter(data["integrated_int_DAPI_norm"], data[col], s=1, alpha=0.1)
            for idx, phase in enumerate(phases):
                phase_df = data.loc[data.cell_cycle == phase]
                ax.scatter(
                    phase_df["integrated_int_DAPI_norm"],
                    phase_df[col],
                    s=1,
                    c=colors[idx],
                    alpha=0.1,
                )
                ax.set_xlabel("norm. DNA content (log2)")
                ax.set_ylim([y_min, y_max])
                ax.set_yscale("log")
                if i == len(condition_list)/2:
                    ax.set_ylabel(col)
                if i == 6:
                    ax.set_ylabel(col)
        ax.set_xscale("log", base=2)
        ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        ax.set_xticks([1, 2, 4, 8])
        ax.set_xlim([1, 16])

        ax.grid(visible=False)
    fig.suptitle(f"{title} {cell_line}", x=0.1, size=16, weight='bold')
    fig._suptitle.set_weight('bold')
    if save:
        save_fig(path, f"{title} {cell_line}", fig_extension="png")

def cellcycle_scatterplot(df: pd.DataFrame, conditions: List[str], cell_line: str, col: str = 'intensity_mean_EdU_nucleus_norm',
                   title='CellCycle Scatter Plot', colors=colors,  y_log=True, save=True, path=Path.cwd()) -> None:
    """
    Function to plot cell cycle scatter plot for a given condition and cell line
    :param df: df from cell_cycle_analysis function
    :param conditions: list of conditions to be plotted
    :param cell_line: cell line to be analysed
    :param col: column to plot at the y-axis, default 'intensity_mean_EdU_nucleus_norm'
    :param title: default 'CellCycle Scatter Plot'
    :param colors: option, default colors from prop_cycle
    :param y_log: boolean, default True, if True plots the y-axis in log scale
    :param save: boolean, default True, if True saves the figure in the path provided
    :param path: option, default: Path.cwd(), path to save the figure
    :return: None, saves figure in path
    """
    df = df[df.cell_line == cell_line]
    phases = ["G1", "S", "G2/M", "Polyploid", "Sub-G1"]
    col_number = len(conditions)
    fig, ax = plt.subplots(ncols=col_number, figsize=(16, 3), sharey='all')
    y_max = df[col].quantile(0.99) * 1.3
    y_min = df[col].quantile(0.01) * 0.8
    for i, condition in enumerate(conditions):
        data = df.loc[df.condition == condition]
        ax[i].scatter(data["integrated_int_DAPI_norm"], data[col], s=1, alpha=0.1)
        for idx, phase in enumerate(phases):
            phase_df = data.loc[data.cell_cycle == phase]

            ax[i].scatter(
                phase_df["integrated_int_DAPI_norm"],
                phase_df[col],
                s=1,
                c=colors[idx],
                alpha=0.1,
            )
        ax[i].set_xlabel("norm. DNA content (log2)")
        ax[i].set_xscale("log", base=2)
        ax[i].xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        ax[i].set_xticks([1, 2, 4, 8])
        ax[i].set_xlim([1, 8])
        ax[i].set_ylim([y_min, y_max])
        ax[i].grid(visible=False)
        ax[i].set_title(conditions[i])
        if i == 0:
            ax[i].set_ylabel(col)
        if y_log:
            ax[i].set_yscale("log")
        fig.suptitle(f"{title} {cell_line}", x=0.1, size=16, weight='bold')
    if save:
        save_fig(path, f"{title} {cell_line}", fig_extension="png")


def cellcycle_histplot(df: pd.DataFrame, conditions: List[str], cell_line: str, title='CellCycle Histogram',
              bins=300, save=True, path=Path.cwd()) -> None:
    """
    Function to plot cell cycle histogram for a given condition and cell line
    :param df: df from cell_cycle_analysis function
    :param conditions:  list of conditions to be plotted
    :param cell_line: cell line to be analysed
    :param title: default 'CellCycle Histogram'
    :param save: boolean, default True, if True saves the figure in the path provided
    :param path: option, default: Path.cwd(), path to save the figure
    :return: None, saves figure in path
    """
    df = df[df.cell_line == cell_line]
    col_number = len(conditions)
    fig, ax = plt.subplots(ncols=col_number, figsize=(16, 3), sharey="all")

    for i, condition in enumerate(conditions):
        data = df.loc[df.condition == condition]
        data["integrated_int_DAPI_norm"].plot.hist(bins=bins, ax=ax[i])
        ax[i].set_xlabel("norm. DNA content (log2)")
        ax[i].set_xscale("log", base=2)
        ax[i].xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        ax[i].set_xticks([1, 2, 4, 8])
        ax[i].set_xlim([1, 16])
        ax[i].grid(visible=False)
        ax[i].set_title(condition)
        if i == 0:
            ax[i].set_ylabel("Frequency")
        fig.suptitle(f"{title} {cell_line}", x=0.1, size=16, weight='bold')
    if save:
        save_fig(path, f"{title}_{cell_line}")

# Proportional cell cycle analysis

def prop_pivot(df_prop: pd.DataFrame, conditions) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Function to pivot the cell cycle proportion dataframe and get the mean and std of each cell cycle phase
    This will be the imput to plot the stacked barplots with errorbars.
    :param df_prop: dataframe from cellcycle_prop function
    :param conditions: list of conditions sort the order of data
    :return: dataframe to submit to the barplot function
    """
    cc_phases = ["G1", "S", "G2/M", "Polyploid", "Sub-G1"]
    df_prop1 = df_prop.copy()
    df_prop1["condition"] = pd.Categorical(df_prop1["condition"], categories=conditions, ordered=True)
    df_mean = df_prop1.groupby(["condition", "cell_cycle"])["percent"].mean().sort_index(level="condition").reset_index().pivot_table(columns=["cell_cycle"], index=["condition"])
    df_mean.columns = df_mean.columns.droplevel(0)
    df_mean = df_mean[cc_phases]
    df_std = df_prop1.groupby(["condition", "cell_cycle"])["percent"].std().sort_index(level="condition").reset_index().pivot_table(columns=["cell_cycle"], index=["condition"])
    df_std.columns = df_std.columns.droplevel(0)
    df_std = df_std[cc_phases]
    return df_mean, df_std

def cellcycle_barplot(df: pd.DataFrame, conditions: List[str],
                      title: str = 'CellCycle Summary Barplot', save: bool = True, path: Path = Path.cwd()) -> None:
    """
    Function to plot cell cycle barplot for a given condition and cell line
    :param dfbar: dataframe from cellcycle_prop function
    :param conditions: list of conditions to be plotted
    :param title: default 'CellCycle Summary Barplot'
    :param save: boolean, default True, if True saves the figure in the path provided
    :param path: option, default: Path.cwd(), path to save the figure
    :return: None, saves figure in path
    """
    cell_lines = df.cell_line.unique()
    fig, ax = plt.subplots(ncols=len(cell_lines), figsize=(6*len(cell_lines), 6))
    if len(cell_lines) == 1:
        df_mean, df_std = prop_pivot(df, conditions)
        df_mean.plot(kind="bar", stacked=True, yerr=df_std, width=0.75, ax=ax)
        ax.set_ylim(0, 110)
        ax.set_xlabel('')
        ax.set_xticklabels(conditions, rotation=30, ha='right')
        ax.set_ylabel('% of population')
        ax.legend(['G1', 'S', 'G2/M', 'Polyploid', 'Sub-G1'], title='CellCyclePhase', bbox_to_anchor=(1, 0.5))

    else:
        for number, cell_line in enumerate(cell_lines):
            df1 = df.loc[df.cell_line == cell_line]
            df_mean, df_std = prop_pivot(df1, conditions)
            df_mean.plot(kind="bar", stacked=True, yerr=df_std, width=0.75, ax=ax[number])
            ax[number].set_ylim(0, 110)
            ax[number].set_xticklabels(conditions, rotation=30, ha='right')
            ax[number].set_xlabel(cell_line)
            if number == 0:
                    ax[number].legend(['G1', 'S', 'G2/M', 'Polyploid', 'Sub-G1'], title='CellCyclePhase', bbox_to_anchor=((len(cell_lines))*1.2, 1))
                    ax[number].set_ylabel('% of population')
            else:
                ax[number].legend().remove()
                ax[number].set_ylabel(None)
    fig.suptitle(title, size=16, weight='bold', x=0.05, horizontalalignment='left')
    if save:
        save_fig(path, title, tight_layout=False)


def intensity_plot(df, condition, conditions, title , path):
    fig, ax = plt.subplots(figsize=(4, 3))

    sns.stripplot(x='condition',
                  y= condition,
                  data = df,
                  order = conditions,
                  size= 2,
                  alpha = 0.3,
                  ax=ax)
    sns.boxplot(x='condition',
                  y= condition,
                  data = df,
                  order= conditions,
                  width = 0.5,
                  ax=ax)
    ax.set_xlabel('')
    ax.set_ylim((0, 20000))
    ax.set_xticklabels(conditions, rotation=30, ha='right')
    fig.suptitle(title, x=0.1, size=16, weight='bold')
    fig._suptitle.set_weight('bold')
    save_fig(path, title)

def violin_plot(df, condition, conditions, title , path):
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.violinplot(x='condition',
                  y= condition,
                  data = df,
                  order = conditions,
                  size= 2,
                  alpha = 0.3,
                  ax=ax)
    ax.set_xlabel('')
    ax.set_xticklabels(conditions, rotation=30, ha='right')
    fig.suptitle(title, ha='left', size=16, weight='bold')
    save_fig(path, title)


def cellcycle_boxplot(df, condition, conditions, colors, y_label, title, path):
    fig, ax = plt.subplots(nrows=3, figsize=(20, 6), sharey=True)
    phases = ['G1', 'S', 'G2/M']
    for row, phase in enumerate(phases):
        sns.boxplot(x='condition', y=condition, data=df[df.cell_cycle==phase], color=colors[row], order=conditions, showfliers=False, ax=ax[row])
        ax[row].set_xlabel('')
        ax[row].set_title(phase)
        if row == 1:
            ax[row].set_ylabel(y_label)
        else:
            ax[row].set_ylabel('')
    fig.suptitle(title, x=0.1, size=16, weight='bold')
    save_fig(path, title)


if __name__ == '__main__':
    plt.style.use("/Users/hh65/matplotlib_style/Style_01.mplstyle")
    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = prop_cycle.by_key()["color"]
    conditions = ['siCtr', 'siCdc27']
    df1 = pd.read_csv('../../data/sample_data_cell_min.csv', index_col=0)
    df_cc = cellcycle_analysis(df1)
    df_cellcycle = cellcycle_prop(df_cc)
    df_bar = barplot_df(df_cellcycle)
    #cellcycle_plot_comb(df_cc, 'intensity_mean_EdU_nucleus_norm', list(df_cc.condition.unique()), colors, 'Comb_Plot', path= pathlib.Path('/Users/hh65/Desktop'))
    #cellcycle_plot(df_cc, 'intensity_mean_EdU_nucleus_norm', list(df_cc.condition.unique()), colors, 'EdU_Plot', path= pathlib.Path('/Users/hh65/Desktop'))
    rel_cellnumber(count_per_cond(df1, 'siCtr'), conditions, colors[:len(df1.condition.unique())], 'test_count', Path('/Users/hh65/Desktop'))
    cellcycle_barplot(df_bar, conditions, 'barplot_test', Path('/Users/hh65/Desktop'))
    intensity_plot(df_cc, 'intensity_mean_EdU_nucleus_norm', conditions,'intensity_plot_test', Path('/Users/hh65/Desktop'))
    # hist_plot(df_cc, df_cc.condition.unique(), 'Hist_Plot', path= Path('/Users/hh65/Desktop'))
