import ifanalysis
import pandas as pd
import pytest
import pandas
@pytest.fixture()
def data():
    return pd.read_csv('../data/sample_data_cell.csv', index_col =0)

def test_agg_multinucleates(data):
    df = ifanalysis.agg_multinucleates(data)
    assert len(df) == 23379
