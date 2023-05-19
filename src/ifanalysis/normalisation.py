import numpy as np
import pandas as pd
from typing import List


# 1) Cell Counting
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

# 2) Cell Cycle Normalisation

norm_colums = ('integrated_int_DAPI', "intensity_mean_EdU_nucleus") # Default columns for cell cycle normalisation
def cellcycle_analysis(df: pd.DataFrame, values: List[str] = norm_colums, cyto: bool = True) -> pd.DataFrame:
    """
    Function to normalise cell cycle data using normalise and assign_ccphase functions for each cell line
    :param df: single cell data from omeroscreen
    :param cyto: True if cytoplasmic data is present
    :return: dataframe with cell cycle and cell cycle detailed columns
    """
    if cyto:
        df_agg = agg_multinucleates(df)
        df_agg_corr = delete_duplicates(df_agg)
    else:
        df_agg_corr = df.copy()
    tempfile = pd.DataFrame()
    for cell_line in df_agg_corr["cell_line"].unique():
        df1 = df_agg_corr.loc[df_agg_corr["cell_line"] == cell_line]
        df_norm = normalise(df1, values)
        df_norm['integrated_int_DAPI_norm'] = df_norm['integrated_int_DAPI_norm'] * 2
        tempfile = pd.concat([tempfile, df_norm])
    return assign_ccphase(data=tempfile)

# Helper Functions for cell cycle normalisation
def agg_multinucleates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Function to aggregate multinucleates by summing up the nucleus area and DAPI intensity
    :param df: single cell data from omeroscreen
    :return: corrected df with aggregated multinucleates
    """
    num_cols = list(df.select_dtypes(include=['float64', 'int64']).columns)
    str_cols = list(df.select_dtypes(include=['object']).columns)
    # define the aggregation functions for each column
    agg_functions = {}
    for col in num_cols:
        if col in ['integrated_int_DAPI', 'area_nucleus']:
            agg_functions[col] = 'sum'
        elif 'max' in col and 'nucleus' in col:
            agg_functions[col] = 'max'
        elif 'min' in col and 'nucleus' in col:
            agg_functions[col] = 'min'
        else:
            agg_functions[col] = 'mean'
    df_agg = df.groupby(str_cols + ['image_id', 'Cyto_ID'], as_index=False).agg(agg_functions)
    df_agg["intensity_mean_EdU_nucleus"] = df_agg["intensity_mean_EdU_nucleus"]/df_agg["intensity_mean_EdU_cyto"]
    return df_agg


def delete_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Function to delete duplicates from the agg_multinucleate dataframe
    :param df: dataframe from agg_multinucleates function
    :return: df with deleted duplicates
    """
    temp_data = pd.DataFrame()
    for image in df["image_id"].unique():
        image_data = df.loc[df.image_id == image].drop_duplicates()
        temp_data = pd.concat([temp_data, image_data])
    return temp_data


def normalise(df: pd.DataFrame, values: list[str]) -> pd.DataFrame:
    """
    Data normalisation function: Identifies the most frequent intensity value and sets it to
    1 by division. For DAPI data this is set to two, to reflect diploid (2N) state of chromosomes
    :param df: dataframe from delete_duplicates function
    :param values:
    :return:
    """
    norm_df = pd.DataFrame()
    for cell_line in df["cell_line"].unique():
        tmp_data = df.copy().loc[(df["cell_line"] == cell_line)]
        tmp_bins = 10000
        for value in values:

            y, x = np.histogram(tmp_data[value], bins=tmp_bins)
            max_value = x[np.where(y == np.max(y))]
            tmp_data[f"{value}_norm"] = tmp_data[value] / max_value[0]
        norm_df = pd.concat([norm_df, tmp_data])
    return norm_df


def assign_ccphase(data: pd.DataFrame) -> pd.DataFrame:
    """
    Assigns a cell cycle phase to each cell based on normalised EdU and DAPI intensities.
    :param data: dataframe from normalise function
    :return: dataframe with cell cycle assignment
    (col: cellcycle (Sub-G1, G1, S, G2/M Polyploid
    and col: cellcycle_detailed with Early S/Late S and Polyploid (non-replicating)
    Polyploid (replicating))
    """
    data["cell_cycle_detailed"] = data.apply(thresholding, axis=1)
    data['cell_cycle'] = data["cell_cycle_detailed"]
    data['cell_cycle'] = data['cell_cycle'].replace(['Early S', 'Late S'], 'S')
    data['cell_cycle'] = data['cell_cycle'].replace(["Polyploid (non-replicating)", "Polyploid (replicating)"],
                                                    'Polyploid')
    return data


def thresholding(data: pd.DataFrame, DAPI_col: str = 'integrated_int_DAPI_norm',
                 EdU_col="intensity_mean_EdU_nucleus_norm") -> str:
    """
    Function to assign cell cycle phase based on thesholds of normalised EdU and DAPI intensities
    :param data: data from assign_ccphase function
    :param DAPI_col: default 'integrated_int_DAPI_norm'
    :param EdU_col: default 'intensity_mean_EdU_nucleus_norm'
    :return: string indicating cell cycle phase
    """
    if data[DAPI_col] <= 1.5:
        return "Sub-G1"

    elif 1.5 < data[DAPI_col] < 3 and data[EdU_col] < 1.5:
        return "G1"

    elif 3 <= data[DAPI_col] < 5.25 and data[EdU_col] < 1.5:
        return "G2/M"

    elif 1.5 < data[DAPI_col] < 3 and data[EdU_col] > 1.5:
        return "Early S"

    elif 3 <= data[DAPI_col] < 5.25 and data[EdU_col] > 1.5:
        return "Late S"

    elif data[DAPI_col] >= 5.25 and data[EdU_col] < 1.5:
        return "Polyploid (non-replicating)"

    elif data[DAPI_col] >= 5.25 and data[EdU_col] > 1.5:
        return "Polyploid (replicating)"

    else:
        return "Unassigned"



# 3 Cell Cycle Proportion Analysis


def cellcycle_prop(df_norm: pd.DataFrame, cell_cycle: str = 'cell_cycle') -> pd.DataFrame:
    """
    Function to calculate the proportion of cells in each cell cycle phase
    :param df_norm: dataframe from assign_ccphase function
    :param cell_cycle: choose column cell_cycle or cell_cycle_detailed, default 'cell_cycle'
    :return: grouped dataframe with cell cycle proportions
    """
    df_ccphase = (
        df_norm.groupby(["plate_id", "well", "cell_line", "condition", cell_cycle])[
            "experiment"
        ].count()
        / df_norm.groupby(["plate_id", "well", "cell_line", "condition"])["experiment"].count()
        * 100
    )
    return df_ccphase.reset_index().rename(columns={"experiment": "percent"})



if __name__ == '__main__':
    df1 = pd.read_csv('../../data/sample_data_cell_min.csv', index_col=0)
    dfagg = agg_multinucleates(df1)
    print(len(dfagg))


    df_cc = cellcycle_analysis(df1)
    print(df_cc.head())

    df_count = count_per_cond(df_cc, 'siCtr')
    print(df_count)

    df_phase = cellcycle_prop(df_cc)
    print(df_phase)
    df_bar = barplot_df(df_phase)
    print(df_bar)


