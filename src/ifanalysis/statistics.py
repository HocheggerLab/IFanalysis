import pandas as pd
from scipy.stats import ttest_ind
from typing import List



def t_test_cc(df: pd.DataFrame, cond: str, cc_phase: str, ctr: str) -> float:
    """
    Function to perform t-test between two conditions for each cell cycle phase
    :param df: dataframe from cellcycle_prop function
    :param cond: condition to compare to control
    :param cc_phase: cell cycle phase to compare
    :param ctr: control condition, default 'SCR'
    :return: p-value
    """
    sample1 = df.loc[(df.cell_cycle == cc_phase) & (df.condition == ctr), 'percent'].to_list()
    sample2 = df.loc[(df.cell_cycle == cc_phase) & (df.condition == cond), 'percent'].to_list()
    t, p = ttest_ind(sample1, sample2)
    return p


# The we loop through each cell cycle ophase using another function
def t_test_perphase(df: pd.DataFrame, cond: str, ctr: str) -> List[float]:
    """
    Function to loop through each cell cycle phase and perform t-test
    :param df: dataframe from cellcycle_prop function
    :param cond: condition to compare to control
    :return:
    """
    return [t_test_cc(df, cond, cc_phase, ctr) for cc_phase in df.cell_cycle.unique()]


def t_test_cc_percond(df, ctr: str ='SCR'):
    """
    Function to loop through each condition and perform t-test for each cell cycle phase
    :param df: dataframe from cellcycle_prop function
    :return: dataframe with mean, std and p-value for each cell cycle phase
    """
    stats_sum = df.groupby(['condition', 'cell_cycle'])['percent'].agg(['mean', 'std']).reset_index()
    p_value_list = []
    for cond in stats_sum.condition.unique():
        p = t_test_perphase(df, cond, ctr)
        p_value_list.append(p)
        # flatten list of lists
    final_pval_list = [item for sublist in p_value_list for item in sublist]
    stats_sum['p_value'] = final_pval_list
    return stats_sum
