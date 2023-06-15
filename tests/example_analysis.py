import matplotlib.pyplot as plt
import pandas as pd
from ifanalysis import normalisation, figures
import pathlib
# Matplotlib Style and Colors
plt.style.use("/Users/hh65/Documents/Current_Coding/hhlab-ifanalysis/tests/matplotlib_style/Style_01.mplstyle")
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]

data_path = pathlib.Path('/Users/hh65/Desktop/OmeroScreen_test/single_cell_data/OmeroScreen_test_final_data.csv')
conditions = ['siCtr', 'siCdc27']
control = 'siCtr'
def default_if_analysis():

    output_path = data_path.parent / 'analysis_files'
    output_path.mkdir(exist_ok=True)

    df = pd.read_csv(data_path, index_col=0)
    df_counts = normalisation.count_per_cond(df, control)
    df_counts.to_csv(output_path / 'counts.csv')
    df_cellcycle = normalisation.cellcycle_analysis(df, H3=True)
    figures.rel_cellnumber(df_counts, conditions, 'U2OS', path=output_path)
    figures.cellcycleplot_comb(df_cellcycle, conditions, 'U2OS', bins=200, path=output_path)
    df_prop = normalisation.cellcycle_prop(df_cellcycle)
    df_prop.to_csv(output_path / 'cellcycle_prop.csv')
    figures.cellcycle_barplot(df_prop, conditions, H3=True, path=output_path)


if __name__ == '__main__':
    default_if_analysis()