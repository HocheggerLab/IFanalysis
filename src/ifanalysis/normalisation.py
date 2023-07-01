import numpy as np
import pandas as pd




norm_colums = ('integrated_int_DAPI', "intensity_mean_EdU_nucleus") # Default columns for cell cycle normalisation
def cellcycle_analysis(df: pd.DataFrame, H3: bool =False, cyto: bool = True) -> pd.DataFrame:
    """
    Function to normalise cell cycle data using normalise and assign_ccphase functions for each cell line
    :param df: single cell data from omeroscreen
    :param cyto: True if cytoplasmic data is present
    :return: dataframe with cell cycle and cell cycle detailed columns
    """
    df1 = df.copy()
    if H3:
        values = ['integrated_int_DAPI', "intensity_mean_EdU_nucleus", "intensity_mean_H3P_nucleus"]
        df1['intensity_mean_H3P_nucleus'] = df1['intensity_mean_H3P_nucleus'] - df1['intensity_min_H3P_nucleus'] + 1
    else:
        values = ['integrated_int_DAPI', "intensity_mean_EdU_nucleus"]
    df1['intensity_mean_EdU_nucleus'] = df1['intensity_mean_EdU_nucleus'] - df1['intensity_min_EdU_nucleus'] + 1
    if cyto:
        df_agg = agg_multinucleates(df1)
        df_agg_corr = delete_duplicates(df_agg)
    else:
        df_agg_corr = df1.copy()
    tempfile = pd.DataFrame()
    for cell_line in df_agg_corr["cell_line"].unique():
        df1 = df_agg_corr.loc[df_agg_corr["cell_line"] == cell_line]
        df_norm = normalise(df1, values)
        df_norm['integrated_int_DAPI_norm'] = df_norm['integrated_int_DAPI_norm'] * 2
        tempfile = pd.concat([tempfile, df_norm])
    return assign_ccphase(data=tempfile, H3=H3)

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
    return df.groupby(str_cols + ['image_id', 'Cyto_ID'], as_index=False).agg(
        agg_functions
    )


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


def assign_ccphase(data: pd.DataFrame, H3) -> pd.DataFrame:
    """
    Assigns a cell cycle phase to each cell based on normalised EdU and DAPI intensities.
    :param data: dataframe from normalise function
    :return: dataframe with cell cycle assignment
    (col: cellcycle (Sub-G1, G1, S, G2/M Polyploid
    and col: cellcycle_detailed with Early S/Late S and Polyploid (non-replicating)
    Polyploid (replicating))
    """
    if H3:
        data["cell_cycle_detailed"] = data.apply(thresholdingH3, axis=1)
    else:
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

    elif 1.5 < data[DAPI_col] < 3 and data[EdU_col] < 3:
        return "G1"

    elif 3 <= data[DAPI_col] < 5.5 and data[EdU_col] < 3:
        return "G2/M"

    elif 1.5 < data[DAPI_col] < 3 and data[EdU_col] > 3:
        return "Early S"

    elif 3 <= data[DAPI_col] < 5.5 and data[EdU_col] > 3:
        return "Late S"

    elif data[DAPI_col] >= 5.5 and data[EdU_col] < 3:
        return "Polyploid (non-replicating)"

    elif data[DAPI_col] >= 5.5 and data[EdU_col] > 3:
        return "Polyploid (replicating)"

    else:
        return "Unassigned"


def thresholdingH3(data: pd.DataFrame, DAPI_col: str = 'integrated_int_DAPI_norm',
                   EdU_col="intensity_mean_EdU_nucleus_norm", H3P_col="intensity_mean_H3P_nucleus_norm") -> str:
    """
    Function to assign cell cycle phase based on thresholds of normalised EdU, DAPI and H3P intensities
    :param data: data from assign_ccphase function
    :param DAPI_col: default 'integrated_int_DAPI_norm'
    :param EdU_col: default 'intensity_mean_EdU_nucleus_norm'
    :param H3P_col: default 'intensity_mean_H3P_nucleus_norm'
    :return: string indicating cell cycle phase
    """
    if data[DAPI_col] <= 1.5:
        return "Sub-G1"

    elif 1.5 < data[DAPI_col] < 3 and data[EdU_col] < 3:
        return "G1"

    elif 3 <= data[DAPI_col] < 5.5 and data[EdU_col] < 3 and data[H3P_col] < 5:
        return "G2"

    elif 3 <= data[DAPI_col] < 5.5 and data[EdU_col] < 3 and data[H3P_col] > 5:
        return "M"

    elif 1.5 < data[DAPI_col] < 3 and data[EdU_col] > 3:
        return "Early S"

    elif 3 <= data[DAPI_col] < 5.5 and data[EdU_col] > 3:
        return "Late S"

    elif data[DAPI_col] >= 5.5 and data[EdU_col] < 3:
        return "Polyploid (non-replicating)"

    elif data[DAPI_col] >= 5.5 and data[EdU_col] > 3:
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



