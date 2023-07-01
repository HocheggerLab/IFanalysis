import pytest
import pandas as pd
import seaborn as sns
from pathlib import Path
from ifanalysis.counts import (
    count_per_cond,
    normalise_count,
    rel_cellnumber,
    abs_cellnumber,
    count_plots,
    count_cells,
)


@pytest.fixture
def sample_df():
    data = {
        'cell_line': ['X', 'X', 'X', 'X', 'Y', 'Y', 'Y', 'Y'],
        'condition': ['A', 'A', 'B', 'B', 'A', 'A', 'B', 'B', ],
        'well': ['W1', 'W2', 'W3', 'W4', 'W5', 'W6', 'W7', 'W8'],
        'plate_id': [1, 1, 1, 1, 1, 1, 1, 1],
        'experiment': [1, 1, 1, 1, 1, 1, 1, 1],
        'abs cell count': [100, 100, 50, 50, 70, 70, 50, 50]
    }
    return pd.DataFrame(data)

@pytest.fixture
def real_df():
    return pd.read_csv("../data/Test_Data_long.csv")

@pytest.mark.parametrize("ctr_cond, expected_mean",
                          [("A", 0.803571), ("B", 1.35)],
                          )
def test_normalise_count_mean(sample_df, ctr_cond, expected_mean):
    mean = normalise_count(sample_df, ctr_cond)[0].mean()
    assert round(mean, 6) == expected_mean
def test_count_per_cond_len(real_df):
    assert len(count_per_cond(real_df, "SCR")) == 96

@pytest.mark.parametrize("column, expected_mean",
                         [("abs cell count", 4572.71875), ("norm_count", 0.67905)],
                         )
def test_count_per_cond_mean(real_df, column, expected_mean):
    mean = count_per_cond(real_df, "SCR")[column].mean()
    assert round(mean, 5) == expected_mean

#
def test_rel_cellnumber(sample_df):
    conditions = ['A', 'B']
    cell_line = 'X'
    title = 'Relative Cell Number'
    save = False
    ax = None
    processed_df = count_per_cond(sample_df[sample_df.cell_line == cell_line], 'A')
    expected_result = sns.barplot(
        data= processed_df,
        x="condition",
        y="norm_count",
        errorbar="sd",
        order=conditions,
        ax=ax,
    )

    result = rel_cellnumber(processed_df, conditions, cell_line, title=title, save=save, ax=ax)

    assert isinstance(result, type(expected_result))
    assert result.get_xlabel() == expected_result.get_xlabel()
    assert result.get_ylabel() == expected_result.get_ylabel()
    assert result.get_title() == expected_result.get_title()

def test_rel_cellnumber_save(sample_df, tmp_path):
    conditions = ['A', 'B']
    cell_line = 'X'
    title = 'Relative Cell Number'
    processed_df = count_per_cond(sample_df[sample_df.cell_line == cell_line], 'A')
    rel_cellnumber(processed_df, conditions, cell_line, title=title, save=True, path=tmp_path)
    result_path= tmp_path / 'Relative Cell Number X.pdf'
    assert result_path.exists()

#
# def test_abs_cellnumber(test_data_frame):
#     conditions = ['A', 'B']
#     cell_line = 'X'
#     title = 'Absolute Cell Number'
#     save = False
#     ax = None
#
#     expected_result = sns.barplot(
#         data=test_data_frame[test_data_frame.cell_line == cell_line],
#         x="condition",
#         y="abs cell count",
#         errorbar="sd",
#         order=conditions,
#         ax=ax,
#     )
#
#     result = abs_cellnumber(test_data_frame, conditions, cell_line, title=title, save=save, ax=ax)
#
#     assert isinstance(result, type(expected_result))
#     assert result.get_xlabel() == expected_result.get_xlabel()
#     assert result.get_ylabel() == expected_result.get_ylabel()
#     assert result.get_title() == expected_result.get_title()
#
#
# def test_count_plots(test_data_frame, tmp_path):
#     conditions = ['A', 'B']
#     title = 'Cell Count Analysis'
#     save = True
#     path = tmp_path / 'test_plot.png'
#
#     count_plots(test_data_frame, conditions, title, save=save, path=path)
#
#     assert path.exists()
#
#
# def test_count_cells(test_data_frame, tmp_path):
#     conditions = ['A', 'B']
#     ctr_cond = 'A'
#     title = 'Cell Count Analysis'
#     save = True
#     path = tmp_path / 'test_plot.png'
#
#     count_cells(test_data_frame, conditions, ctr_cond, title=title, save=save, path=path)
#
#     assert path.exists()
