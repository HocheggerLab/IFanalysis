import numpy as np
import pandas as pd



def agg_multinucleates(df):
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
    #df_agg["intensity_mean_EdU_nucleus"] = df_agg["intensity_mean_EdU_nucleus"]/df_agg["intensity_mean_EdU_cyto"]
    return df_agg


def delete_duplicates(df):
    """delete duplicate cell_IDs for each image"""
    temp_data = pd.DataFrame()
    for image in df["image_id"].unique():
        image_data = df.loc[df.image_id == image].drop_duplicates()
        temp_data = pd.concat([temp_data, image_data])
    return temp_data


def normalise(df: pd.DataFrame, values: list[str]) -> pd.DataFrame:
    """
    Data normalisation function: Identifies the most frequent intensity value and sets it to
    1 by division. For DAPI data this is set to two, to reflect diploid (2N) state of chromosomes
    """
    tmp_output = pd.DataFrame()
    for cell_line in df["cell_line"].unique():
        tmp_data = df.copy().loc[(df["cell_line"] == cell_line)]
        tmp_bins = 10000
        for value in values:

            y, x = np.histogram(tmp_data[value], bins=tmp_bins)
            max_value = x[np.where(y == np.max(y))]
            tmp_data[f"{value}_norm"] = tmp_data[value] / max_value[0]
        tmp_output = pd.concat([tmp_output, tmp_data])
    return tmp_output


def assign_ccphase(data: pd.DataFrame) -> pd.DataFrame:
    """
    Assigns a cell cycle phase to each cell based on normalised EdU and DAPI intensities.
    :param data:
    :return: dataframe with cell cycle assignment (col: cellcycle and col: cellcycle_detailed)
    """
    data["cell_cycle_detailed"] = data.apply(thresholding, axis=1)
    data['cell_cycle'] = data["cell_cycle_detailed"]
    data['cell_cycle'] = data['cell_cycle'].replace(['Early S', 'Late S'], 'S')
    data['cell_cycle'] = data['cell_cycle'].replace(["Polyploid (non-replicating)", "Polyploid (replicating)"],
                                                    'Polyploid')
    return data


def thresholding(
        data: pd.DataFrame, DAPI_col: str = 'integrated_int_DAPI_norm', EdU_col="intensity_mean_EdU_nucleus_norm"
):
    if data[DAPI_col] <= 1.5:
        return "Sub-G1"

    elif 1.5 < data[DAPI_col] < 3 and data[EdU_col] < 2:
        return "G1"

    elif 3 <= data[DAPI_col] < 5.25 and data[EdU_col] < 2:
        return "G2/M"

    elif 1.5 < data[DAPI_col] < 3 and data[EdU_col] > 2:
        return "Early S"

    elif 3 <= data[DAPI_col] < 5.25 and data[EdU_col] > 2:
        return "Late S"

    elif data[DAPI_col] >= 5.25 and data[EdU_col] < 2:
        return "Polyploid (non-replicating)"

    elif data[DAPI_col] >= 5.25 and data[EdU_col] > 2:
        return "Polyploid (replicating)"

    else:
        return "Unassigned"


def cellcycle_analysis(df):
    df_agg = agg_multinucleates(df)
    df_norm = normalise(
        df=df_agg,
        values=['integrated_int_DAPI', "intensity_mean_EdU_nucleus"])
    df_norm['integrated_int_DAPI_norm'] = df_norm['integrated_int_DAPI_norm'] * 2
    return assign_ccphase(data=df_norm)

def count_percond(df, ctr_cond):
    df_count = (
        df.groupby(["cell_line", "condition", "well", "plate_id"])["experiment"]
        .count()
        .reset_index()
        .rename(columns={"experiment": "abs cell count"})
    )
    df_count["norm_count"] = normalise_cell_line(df_count, ctr_cond)
    return df_count

def normalise_cell_line(df, ctr_cond):
    norm_count = pd.DataFrame()
    for cell_line in df["cell_line"].unique():
        norm_value = df.loc[
            (df["cell_line"] == cell_line) & (df["condition"] == ctr_cond),
            "abs cell count",
        ].mean()
        rel_cellcount = (
                df.loc[df["cell_line"] == cell_line, "abs cell count"] / norm_value
        )
        return pd.concat([norm_count, rel_cellcount])

def rel_cellcycle_phase(df, cell_cycle = 'cell_cycle'):
    df_ccphase = (
        df.groupby(["plate_id", "well", "cell_line", "condition", cell_cycle])[
            "experiment"
        ].count()
        / df.groupby(["plate_id", "well", "cell_line", "condition"])["experiment"].count()
        * 100
    )
    return df_ccphase.reset_index().rename(columns={"experiment": "percent"})

def barplot_df(df_cellcycle):
    dfbar = (
        df_cellcycle.groupby(["cell_line", "condition", "cell_cycle"])["percent"]
        .mean()
        .reset_index()
        .pivot_table(columns=["cell_cycle"], index=["condition", "cell_line"])
        .reset_index()

    )
    col_names = ['condition', 'cell_line', 'Sub-G1', 'G1', 'S', 'G2/M', 'Polypolid']

    dfbar.columns = dfbar.columns.map(lambda x: ''.join(x))

    dfbar = dfbar.reindex(columns=['condition', 'cell_line', 'percentSub-G1', 'percentG1', 'percentS', 'percentG2/M', 'percentPolyploid'])
    dfbar = dfbar.rename(columns=dict(zip(dfbar.columns, col_names)))
    return dfbar



if __name__ == '__main__':
    df1 = pd.read_csv('../../data/sample_data_cell_min.csv', index_col=0)
    dfagg = agg_multinucleates(df1)
    print(len(dfagg))


    df_cc = cellcycle_analysis(df1)
    print(df_cc.head())

    df_count = count_percond(df_cc, 'siCtr')
    print(df_count)

    df_phase = rel_cellcycle_phase(df_cc)
    print(df_phase)
    df_bar = barplot_df(df_phase)
    print(df_bar)


