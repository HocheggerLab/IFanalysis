import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
import seaborn as sns
import pathlib
from ifanalysis.normalisation import *


def save_fig(
    path: pathlib, fig_id: str, tight_layout=True, fig_extension="pdf", resolution=300
) -> None:
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
    # print("Saving figure", fig_id)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(dest, format=fig_extension, dpi=resolution)

def cellcycle_plot(df, col, conditions, colors, title, y_log=True, save=True, path=pathlib.Path.cwd()):
    phases = ["G1", "S", "G2/M"]
    col_number = len(conditions)

    fig, ax = plt.subplots(ncols=col_number, figsize=(16, 3), sharey='all')

    for i, condition in enumerate(conditions):
        data = df.loc[df.condition == condition]
        y_max = data[col].quantile(0.99) * 1.2
        y_min = data[col].quantile(0.01) * 0.8

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
        fig.suptitle(title, x=0.1, size=16, weight='bold')
    if save:
        save_fig(path, title)

def cellcycle_plot_comb(df, conditions, colors, title, col='intensity_mean_EdU_nucleus_norm', save=True, path=pathlib.Path.cwd()):
    phases = ["G1", "S", "G2/M"]
    col_number = len(conditions)

    condition_list = conditions * 2
    fig = plt.figure(figsize=(3*len(conditions), 5))
    # define the grid layout with different height ratios
    gs = GridSpec(2, col_number, height_ratios=[1, 3])
    ax_list = [(i,j) for i in range(2) for j in range(col_number)]
    for i, pos in enumerate(ax_list):
        data = df[df.condition == condition_list[i]]
        y_max = data[col].quantile(0.99) * 1.3
        y_min = data[col].quantile(0.01) * 0.8
        ax = fig.add_subplot(gs[pos[0], pos[1]])
        if i < len(condition_list)/2:
            data['integrated_int_DAPI_norm'].plot.hist(bins=500, density=True, color=colors[-1], ax=ax)
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
    fig.suptitle(title, x=0.1, size=16, weight='bold')
    fig._suptitle.set_weight('bold')

    if save:
        save_fig(path, title)

def rel_cellnumber(df_count, conditions, colors, title, path):
    fig, ax = plt.subplots(figsize=(15, 5))
    sns.barplot(
        data=df_count,
        x="condition",
        y="norm_count",
        errorbar="sd",
        order=conditions,
        palette=colors,
        ax=ax,
    )
    ax.set_title(title)
    save_fig(path, title, fig_extension="pdf")


def hist_plot(df, conditions, title, save=True, path=pathlib.Path.cwd()):
    col_number = len(conditions)
    fig, ax = plt.subplots(ncols=col_number, figsize=(16, 3), sharey="all")

    for i, condition in enumerate(conditions):
        data = df.loc[df.condition == condition]
        data["integrated_int_DAPI_norm"].plot.hist(bins=500, ax=ax[i])
        ax[i].set_xlabel("norm. DNA content (log2)")
        ax[i].set_xscale("log", base=2)
        ax[i].xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        ax[i].set_xticks([1, 2, 4, 8])
        ax[i].set_xlim([1, 16])
        ax[i].grid(visible=False)
        ax[i].set_title(condition)
        if i == 0:
            ax[i].set_ylabel("Frequency")
    if save:
        save_fig(path, title)

def cellcycle_barplot(dfbar, conditions, title, path):
    cell_lines = dfbar.cell_line.unique()
    dfbar = dfbar.set_index('condition')
    dfbar = dfbar.reindex(conditions)
    dfbar = dfbar.reset_index()
    fig, ax = plt.subplots(ncols=len(cell_lines), figsize=(6*len(cell_lines), 6), sharey="all")
    if len(cell_lines) == 1:
        dfbar.plot(x="condition", kind="bar", stacked=True, width=0.75, ax=ax)
        ax.set_xlabel('')
        ax.set_xticklabels(conditions, rotation=30, ha='right')
        ax.set_ylabel('% of population')
        ax.legend(['Sub-G1', 'G1', 'S', 'G2/M', 'Polyploid'], title='CellCyclePhase', bbox_to_anchor=(1.3, 0.8))

    else:
        for number, cell_line in enumerate(cell_lines):
            dfbar.loc[dfbar.cell_line == cell_line].plot(x="condition", kind="bar", stacked=True, width=0.75, ax=ax[number])
            ax.set_xticklabels(conditions, rotation=30, ha='right')
            ax[number].set_xlabel(cell_line)
            if number == 0:
                    ax[number].legend(['Sub-G1', 'G1', 'S', 'G2/M', 'Polyploid'], title='CellCyclePhase', bbox_to_anchor=((len(cell_lines))+1.2, 1.2))
                    ax[number].set_ylabel('% of population')
            else:
                ax[number].legend().remove()
                ax[number].set_ylabel(None)
    fig.suptitle(title, x=0.1, size=16, weight='bold')
    fig._suptitle.set_weight('bold')
    save_fig(path, title, fig_extension="pdf")


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
    df_cellcycle = rel_cellcycle_phase(df_cc)
    df_bar = barplot_df(df_cellcycle)
    #cellcycle_plot_comb(df_cc, 'intensity_mean_EdU_nucleus_norm', list(df_cc.condition.unique()), colors, 'Comb_Plot', path= pathlib.Path('/Users/hh65/Desktop'))
    #cellcycle_plot(df_cc, 'intensity_mean_EdU_nucleus_norm', list(df_cc.condition.unique()), colors, 'EdU_Plot', path= pathlib.Path('/Users/hh65/Desktop'))
    rel_cellnumber(count_percond(df1, 'siCtr'), conditions, colors[:len(df1.condition.unique())], 'test_count', pathlib.Path('/Users/hh65/Desktop'))
    cellcycle_barplot(df_bar, conditions, 'barplot_test', pathlib.Path('/Users/hh65/Desktop'))
    intensity_plot(df_cc, 'intensity_mean_EdU_nucleus_norm', conditions,'intensity_plot_test', pathlib.Path('/Users/hh65/Desktop'))
    # hist_plot(df_cc, df_cc.condition.unique(), 'Hist_Plot', path= pathlib.Path('/Users/hh65/Desktop'))
