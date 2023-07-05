import matplotlib.pyplot as plt
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

def standard_cellcycleplots(df: pd.DataFrame, conditions: List[str], H3=False, file_path=path) -> None:
    df_cc = cellcycle_analysis(df, H3=H3)
    for cell_line in df_cc.cell_line.unique():
        cellcycleplot_comb(df_cc, conditions, cell_line, H3=H3, path=file_path)
    df_prop = cellcycle_prop(df_cc)
    cellcycle_barplot(df_prop, conditions, H3=H3, path=file_path)

def cellcycleplot_comb(df_cc: pd.DataFrame, conditions: List[str], cell_line: str, H3: bool = False, bins=1000,
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
    if H3:
        phases = ["G1", "S", "G2", "M", "Polyploid", "Sub-G1"]
    else:
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
            # Histogram Plot on the top row
            data['integrated_int_DAPI_norm'].plot.hist(bins=bins, color=colors[-1], ax=ax)
            ax.set_title(condition_list[i])
            ax.set_xlabel(None)
            if i == 0:
                ax.set_ylabel("Frequency")
            else:
                ax.set_ylabel('')
        else:
            # Scatter Plot on the bottom row
            ax.scatter(data["integrated_int_DAPI_norm"], data[col], s=1, alpha=0.1)
            # Create an empty list to store the scatter plot objects
            scatter_plots = []
            for idx, phase in enumerate(phases):
                phase_df = data.loc[data.cell_cycle == phase]
                scatter = ax.scatter(
                    phase_df["integrated_int_DAPI_norm"],
                    phase_df[col],
                    s=1,
                    c=colors[idx],
                    alpha=0.1,
                    label=phase  # Add a label for the legend
                )
                ax.set_xlabel("norm. DNA content (log2)")
                ax.set_ylim([y_min, y_max])
                ax.set_yscale("log")
                if i == len(condition_list)/2:
                    ax.set_ylabel(col)
                if i == 6:
                    ax.set_ylabel(col)
                scatter_plots.append(scatter)
            # Add legend using the scatter plots and their labels
            if i == len(condition_list)/2:
                legend = ax.legend(handles=scatter_plots, loc="upper right")

                # Increase the alpha (transparency) of the legend entries
                for lh in legend.legendHandles:
                    lh.set_alpha(0.5)
                    lh.set_sizes([10])

                # Adjust the appearance of the legend
                frame = legend.get_frame()
                frame.set_facecolor('white')
                frame.set_edgecolor('black')


        ax.set_xscale("log", base=2)
        ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        ax.set_xticks([1, 2, 4, 8])
        ax.set_xlim([1, 16])

        ax.grid(visible=False)   
        
    fig.suptitle(title, size=16, weight='bold', x=0.05, horizontalalignment='left')
    fig._suptitle.set_weight('bold')
    if save:
        save_fig(path, f"{title} {cell_line}", fig_extension="pdf")

def cellcycle_scatterplot(df: pd.DataFrame, conditions: List[str], cell_line: str, H3: bool = False, col: str = 'intensity_mean_EdU_nucleus_norm',
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
    if H3:
        phases = ["G1", "S", "G2", "M", "Polyploid", "Sub-G1"]
    else:
        phases = ["G1", "S", "G2/M", "Polyploid", "Sub-G1"]
    col_number = len(conditions)

    fig, ax = plt.subplots(ncols=col_number, figsize=(3*len(conditions), 4), sharey='all')
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
        ax[i].set_xlim([1, 16])
        ax[i].set_ylim([y_min, y_max])
        ax[i].grid(visible=False)
        ax[i].set_title(conditions[i])
        if i == 0:
            ax[i].set_ylabel(col)
        if y_log:
            ax[i].set_yscale("log")
        fig.suptitle(f"{title} {cell_line}", size=16, weight='bold', x=0.05, horizontalalignment='left')
    if save:
        save_fig(path, f"{title} {cell_line}", fig_extension="pdf")


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

def prop_pivot(df_prop: pd.DataFrame, conditions, H3=False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Function to pivot the cell cycle proportion dataframe and get the mean and std of each cell cycle phase
    This will be the imput to plot the stacked barplots with errorbars.
    :param df_prop: dataframe from cellcycle_prop function
    :param conditions: list of conditions sort the order of data
    :param H3: boolean, default False, if True the function will use M phase instead of G2/M based on H3 staining
    :return: dataframe to submit to the barplot function
    """
    if H3:
        cc_phases = ["G1", "S", "G2", "M", "Polyploid", "Sub-G1"]
    else:   
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

def cellcycle_barplot(df: pd.DataFrame, conditions: List[str], H3 = False,
                      title: str = 'Cell Cycle Summary Barplot', save: bool = True, path: Path = Path.cwd()) -> None:
    """
    Function to plot cell cycle barplot for a given condition and cell line
    :param df: dataframe from cellcycle_analysis function
    :param conditions: list of conditions to be plotted
    :param title: default 'CellCycle Summary Barplot'
    :param save: boolean, default True, if True saves the figure in the path provided
    :param path: option, default: Path.cwd(), path to save the figure
    :return: None, saves figure in path
    """
    cell_lines = df.cell_line.unique()
    fig, ax = plt.subplots(ncols=len(cell_lines), figsize=(5*len(cell_lines), 5))
    if len(cell_lines) == 1:
        df_mean, df_std = prop_pivot(df, conditions, H3)
        df_mean.plot(kind="bar", stacked=True, yerr=df_std, width=0.75, ax=ax)
        ax.set_ylim(0, 110)
        ax.set_xlabel('')
        ax.set_xticklabels(conditions, rotation=30, ha='right')
        ax.set_ylabel('% of population')
        if H3:
            ax.legend(['G1', 'S', 'G2', 'M', 'Polyploid', 'Sub-G1'], title='CellCyclePhase', bbox_to_anchor=(1, 0.5))
        else:
            ax.legend(['G1', 'S', 'G2/M', 'Polyploid', 'Sub-G1'], title='CellCyclePhase', bbox_to_anchor=(1, 0.5))

    else:
        for number, cell_line in enumerate(cell_lines):
            df1 = df.loc[df.cell_line == cell_line]
            df_mean, df_std = prop_pivot(df1, conditions)
            df_mean.plot(kind="bar", stacked=True, yerr=df_std, width=0.75, ax=ax[number])
            ax[number].set_ylim(0, 110)
            ax[number].set_xticklabels(conditions, rotation=30, ha='right')
            ax[number].set_xlabel(cell_line)
            if number == 0 and H3:
                    ax[number].legend(['G1', 'S', 'G2', 'M', 'Polyploid', 'Sub-G1'], title='CellCyclePhase', bbox_to_anchor=((len(cell_lines))*1.2, 1))
                    ax[number].set_ylabel('% of population')
            elif number == 0:
                    ax[number].legend(['G1', 'S', 'G2/M', 'Polyploid', 'Sub-G1'], title='CellCyclePhase', bbox_to_anchor=((len(cell_lines))*1.2, 1))
                    ax[number].set_ylabel('% of population')
            else:
                ax[number].legend().remove()
                ax[number].set_ylabel(None)
    fig.suptitle(title, size=16, weight='bold', x=0.05, horizontalalignment='left')
    if save:
        save_fig(path, title, tight_layout=False)



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