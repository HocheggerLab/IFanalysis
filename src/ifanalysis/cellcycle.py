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
module_dir = Path(__file__).parent
style_path = module_dir / 'styles/Style_01.mplstyle'
plt.style.use(str(style_path))
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]


path = Path.cwd()



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
    df_mean = df_prop1.groupby(["condition", "cell_cycle"])["percent"].mean().sort_index(level="condition").reset_index().pivot_table(columns=["cell_cycle"], index=["condition"], observed=False)
    df_mean.columns = df_mean.columns.droplevel(0)
    df_mean = df_mean[cc_phases]
    df_std = df_prop1.groupby(["condition", "cell_cycle"])["percent"].std().sort_index(level="condition").reset_index().pivot_table(columns=["cell_cycle"], index=["condition"], observed=False)
    df_std.columns = df_std.columns.droplevel(0)
    df_std = df_std[cc_phases]
    return df_mean, df_std

def cellcycle_barplot(df: pd.DataFrame, conditions: List[str], title_str, H3 = False, save: bool = True, path: Path = Path.cwd()) -> None:
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
        ax=[ax]
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
    fig.suptitle(title_str, size=16, weight='bold', x=0.05, horizontalalignment='left')
    if save:
        save_fig(path, f"barplot_{title_str}", tight_layout=False)






