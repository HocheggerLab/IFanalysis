import matplotlib.pyplot as plt
import pandas as pd
from ifanalysis import normalisation, cellcycle
import pathlib


data_path = pathlib.Path('../data/Test_Data.csv')
conditions = ['siCtr', 'siCdc27']
control = 'siCtr'
def default_if_analysis():

    output_path = data_path.parent / 'sample_figures'
    output_path.mkdir(exist_ok=True)
    df = pd.read_csv(data_path, index_col=0)
    df_counts = normalisation.count_per_cond(df, control)
    df_counts.to_csv(output_path / 'counts.csv')
    df_cellcycle = normalisation.cellcycle_analysis(df, H3=True)
    figures.rel_cellnumber(df_counts, conditions, 'U2OS', path=output_path)
    figures.cellcycleplot_comb(df_cellcycle, conditions, 'U2OS', H3= True, bins=200, path=output_path)
    df_prop = normalisation.cellcycle_prop(df_cellcycle)
    df_prop.to_csv(output_path / 'cellcycle_prop.csv')
    figures.cellcycle_barplot(df_prop, conditions, H3=True, path=output_path)


if __name__ == '__main__':
    default_if_analysis()