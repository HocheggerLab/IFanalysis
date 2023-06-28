import ifanalysis
import pandas as pd
import pytest
import pandas
@pytest.fixture()
def data():
    return pd.read_csv('../data/Test_Data.csv', index_col =0)

def test_agg_multinucleates_len(data):
    df = ifanalysis.agg_multinucleates(data)
    assert len(df) == 26407
    assert len(df) != 29327

def test_agg_multinucleates_cols(data):
    df = ifanalysis.agg_multinucleates(data)
    assert df.columns[1:4].tolist() == ['well', 'cell_line', 'condition']

@pytest.mark.parametrize(['column', 'result'],
                         [('intensity_mean_DAPI_nucleus', 16214.864537457832),
                          ('intensity_mean_DAPI_cyto', 2490.507574899092),
                          ('integrated_int_DAPI', 2325766.4839166077),
                          ('area_nucleus', 159.35835195213392)
                          ])

def test_agg_multinucleates_mean(data, column, result):
    df = ifanalysis.agg_multinucleates(data)
    assert df[column].mean() == result


