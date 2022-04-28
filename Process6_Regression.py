# Process6_Regression.py

__author__ = "Scott Matsumura"
__last_update__ = "4.15.2022"

import os
import math
import copy
import numpy as np
from scipy import stats
import pandas as pd
from statsmodels.compat import lzip
import statsmodels.api as sm
import statsmodels.stats.api as sms
import matplotlib.pyplot as plt

###################
""" Data options """

LIST_OPTIONS = [
    {'independent_data': ['amino_entropy'],
     'dependent_data': 'divergence_time',
     'label': 'AE_All',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['protein_entropy'],
     'dependent_data': 'divergence_time',
     'label': 'PE_All',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['constructive_entropy'],
     'dependent_data': 'divergence_time',
     'label': 'CE_All',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['amino_entropy'],
     'dependent_data': 'divergence_time',
     'taxonomy_rank': 'superkingdom_name',
     'taxonomy_name': 'Eukaryota',
     'label': 'AE_Eukaryota',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['protein_entropy'],
     'dependent_data': 'divergence_time',
     'taxonomy_rank': 'superkingdom_name',
     'taxonomy_name': 'Eukaryota',
     'label': 'PE_Eukaryota',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['constructive_entropy'],
     'dependent_data': 'divergence_time',
     'taxonomy_rank': 'superkingdom_name',
     'taxonomy_name': 'Eukaryota',
     'label': 'CE_Eukaryota',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['amino_entropy'],
     'dependent_data': 'divergence_time',
     'taxonomy_rank': 'superkingdom_name',
     'taxonomy_name': 'Bacteria',
     'label': 'AE_Bacteria',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['protein_entropy'],
     'dependent_data': 'divergence_time',
     'taxonomy_rank': 'superkingdom_name',
     'taxonomy_name': 'Bacteria',
     'label': 'PE_Bacteria',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['constructive_entropy'],
     'dependent_data': 'divergence_time',
     'taxonomy_rank': 'superkingdom_name',
     'taxonomy_name': 'Bacteria',
     'label': 'CE_Bacteria',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['amino_entropy'],
     'dependent_data': 'divergence_time',
     'taxonomy_rank': 'kingdom_name',
     'taxonomy_name': 'Metazoa',
     'label': 'AE_Metazoa',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['protein_entropy'],
     'dependent_data': 'divergence_time',
     'taxonomy_rank': 'kingdom_name',
     'taxonomy_name': 'Metazoa',
     'label': 'PE_Metazoa',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['constructive_entropy'],
     'dependent_data': 'divergence_time',
     'taxonomy_rank': 'kingdom_name',
     'taxonomy_name': 'Metazoa',
     'label': 'CE_Metazoa',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['amino_entropy'],
     'dependent_data': 'divergence_time',
     'taxonomy_rank': 'kingdom_name',
     'taxonomy_name': 'Viridiplantae',
     'label': 'AE_Viridiplantae',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['protein_entropy'],
     'dependent_data': 'divergence_time',
     'taxonomy_rank': 'kingdom_name',
     'taxonomy_name': 'Viridiplantae',
     'label': 'PE_Viridiplantae',
     'divergence_time': [lambda x: math.log(x, 10), 'log']},
    {'independent_data': ['constructive_entropy'],
     'dependent_data': 'divergence_time',
     'taxonomy_rank': 'kingdom_name',
     'taxonomy_name': 'Viridiplantae',
     'label': 'CE_Viridiplantae',
     'divergence_time': [lambda x: math.log(x, 10), 'log']}]

###################

# File location
PATH_FASTA = '/My Drive/PyCharm/DNAEntropy/fasta/'
DIRECTORY_LOAD = 'zz_CSVdata'
FILE_LOAD = 'entropy_data.csv'

# Select alternative sequences when available
ORTHODB_SEQ = 0

# Display regression plot
DISPLAY_PLOT = True

###################


def apply_regression(path_fasta, directory_load, file_load, orthodb_seq, list_options, display_plot):
    """ Applies a linear regression on the data.

    :param path_fasta: Location of fasta directory from the home directory, which will also be used as save directory.
    :param directory_load: Directory of the location of the data csv files to be used.
    :param file_load: File name of the data csv to be used.
    :param orthodb_seq: OrthoDB sequence number to be used. '0' will use first available sequence.
        '1' will choose the next if available and so on.
    :param list_options: List of selectable dictionary options.
        key: 'independent_data', value: List of key names of the independent data.
        key: 'dependent_data', value: Name for dependent data.
        key: 'taxonomy_rank', value: Taxonomy Rank.
        key: 'taxonomy_name', value: Name within Taxonomy Rank to select out of data.
        key: 'label', value: The will be used as a label for the linear regression and as a file name extension for the
            csv file.
        key: (data name), value: List of two elements. The first is the lambda transformation used. The second is
            the name extension for the label.
    :param display_plot: If True, it will display the linear regression scatter plot.
    """

    # Directory and file names
    dir_csvregression = 'zz_CSVregression'
    file_csvregression = 'entropy_regression'
    file_df_data_div = 'entropy_data'

    path_home = os.getenv('HOME')
    path_save = path_home + path_fasta + '/' + dir_csvregression + '/'

    if not os.path.exists(path_save):
        os.makedirs(path_save)

    path_load_data = path_home + path_fasta + '/' + directory_load + '/' + file_load
    df_data = pd.read_csv(path_load_data)

    # Remove sequences without a divergence time
    df_data_div = df_data.loc[np.isfinite(df_data['divergence_time'])].reset_index(drop=True)

    # Remove duplicate sequences
    list_ncbi_id = list(df_data_div['NCBI_id'])
    list_orthodb_seq = list(df_data_div['OrthoDB_seq'])
    list_duplicates = []
    dict_duplicates_index = {}
    dict_duplicates_seq = {}
    for int_index, (int_ncbi_id, int_orthodb_seq) in enumerate(zip(list_ncbi_id, list_orthodb_seq)):
        if dict_duplicates_index.get(int_ncbi_id) is None:
            dict_duplicates_index[int_ncbi_id] = [int_index]
            dict_duplicates_seq[int_ncbi_id] = [int_orthodb_seq]
        else:
            if len(dict_duplicates_index[int_ncbi_id]) == 1:
                list_duplicates.append(int_ncbi_id)
            dict_duplicates_index[int_ncbi_id].append(int_index)
            dict_duplicates_seq[int_ncbi_id].append(int_orthodb_seq)

    for int_dup in list_duplicates:
        list_dup_index = dict_duplicates_index[int_dup]
        list_dup_seq = dict_duplicates_seq[int_dup]
        list_dup_seq, list_dup_index = (list(t) for t in zip(*sorted(zip(list_dup_seq, list_dup_index))))
        dict_duplicates_index[int_dup] = list_dup_index
        dict_duplicates_seq[int_dup] = list_dup_seq

    dict_alt_counts = {}
    for int_key, list_value in dict_duplicates_index.items():
        int_size_value = len(list_value)
        if dict_alt_counts.get(int_size_value) is None:
            dict_alt_counts[int_size_value] = 1
        else:
            dict_alt_counts[int_size_value] = dict_alt_counts[int_size_value] + 1
    list_keys_counts = list(dict_alt_counts.keys())
    list_keys_counts.sort()

    list_deletions = []
    for int_ncbi_id in list_duplicates:
        list_index = dict_duplicates_index[int_ncbi_id]
        int_size_index = len(list_index)
        if orthodb_seq > int_size_index - 1:
            del list_index[int_size_index - 1]
        else:
            del list_index[orthodb_seq]
        list_deletions = list_deletions + list_index

    df_data_div.drop(df_data_div.index[list_deletions], inplace=True)
    df_data_div.reset_index(drop=True, inplace=True)

    # Count of alternate to concatenate to CSV Stats
    df_data_seq_alt_name = pd.DataFrame({0: 'Counts of alternate', 1: 'species sequences'}, index=[0])
    df_data_seq_alt = pd.DataFrame(dict_alt_counts, columns=list_keys_counts, index=[0])
    df_data_seq_alt = pd.concat([pd.DataFrame({'duplicates': 'counts'}, index=[0]), df_data_seq_alt], axis=1)
    df_data_seq_alt = df_data_seq_alt.T.reset_index().T.reset_index(drop=True)

    # Begin processing options
    for dict_options in list_options:
        # Apply transformations
        list_independent = dict_options['independent_data']
        str_dependent = dict_options['dependent_data']
        str_descriptor_name = dict_options['label']
        str_taxonomy_rank = dict_options.get('taxonomy_rank')
        str_taxonomy_name = dict_options.get('taxonomy_name')
        list_basic_stats_names = [str_dependent] + list_independent

        list_independent_new = []
        str_dependent_new = str_dependent

        if str_taxonomy_rank is not None:
            df_data_process = df_data_div[df_data_div[str_taxonomy_rank] == str_taxonomy_name].copy(deep=True)
        else:
            df_data_process = df_data_div.copy(deep=True)

        # Transformation option for dependent data
        if dict_options.get(str_dependent_new) is not None:
            str_dependent_new = str_dependent + '_' + dict_options[str_dependent][1]
            lambda_transform = dict_options[str_dependent][0]
            list_data_new = list(df_data_process[str_dependent])
            list_data_new_tr = [lambda_transform(x) for x in list_data_new]
            df_data_process[str_dependent_new] = list_data_new_tr
            list_basic_stats_names += [str_dependent_new]

        # Transformation option for independent data
        for str_independent in list_independent:
            if dict_options.get(str_independent) is None:
                list_independent_new.append(str_independent)
            else:
                str_independent_new = str_independent + '_' + dict_options[str_independent][1]
                list_independent_new.append(str_independent_new)
                lambda_transform = dict_options[str_independent][0]
                list_data_new = list(df_data_process[str_independent])
                list_data_new_tr = [lambda_transform(x) for x in list_data_new]
                df_data_process[str_independent_new] = list_data_new_tr
                list_basic_stats_names += [str_independent_new]

        # Apply linear regression
        df_regression = linear_regression(df_data=df_data_process,
                                          list_independent=list_independent_new,
                                          str_dependent=str_dependent_new,
                                          str_descriptor_name=str_descriptor_name,
                                          display_plot=display_plot)

        # Concatenates information to CSV Stats
        # Data: count, minimum, maximum, mean, standard deviation
        list_basic_stats_column_names = ['names',
                                         'count',
                                         'minimum',
                                         'maximum',
                                         'mean',
                                         'standard deviation']
        list_basic_count = []
        list_basic_minimum = []
        list_basic_maximum = []
        list_basic_mean = []
        list_basic_sd = []
        for str_index_name in list_basic_stats_names:
            list_basic_count.append(df_data_process[str_index_name].count())
            list_basic_minimum.append(df_data_process[str_index_name].min())
            list_basic_maximum.append(df_data_process[str_index_name].max())
            list_basic_mean.append(df_data_process[str_index_name].mean())
            list_basic_sd.append(df_data_process[str_index_name].std())
        list_basic_stats_zip = list(zip(list_basic_stats_names,
                                        list_basic_count,
                                        list_basic_minimum,
                                        list_basic_maximum,
                                        list_basic_mean,
                                        list_basic_sd))
        df_basic_stats = pd.DataFrame(list_basic_stats_zip, columns=list_basic_stats_column_names)
        df_basic_stats = df_basic_stats.T.reset_index().T.reset_index(drop=True)
        df_basic_stats_name = pd.DataFrame({0: 'Basic Statistics'}, index=[0])

        # Concatenate CSV Stats to save
        df_empty_series = pd.Series([''])
        df_regression = pd.concat([df_regression,
                                   df_empty_series,
                                   df_empty_series,
                                   df_empty_series,
                                   df_data_seq_alt_name,
                                   df_data_seq_alt,
                                   df_empty_series,
                                   df_basic_stats_name,
                                   df_basic_stats], ignore_index=True, sort=False)

        # Save to CSV
        # CSV Stats
        path_save_regression =\
            path_save + file_csvregression + '_' + str(orthodb_seq) + '_' + str_descriptor_name + '.csv'
        df_regression.to_csv(path_save_regression, header=False, index=False)

        # CSV Data
        path_save_data_div = path_save + file_df_data_div + '_' + str(orthodb_seq) + '_' + str_descriptor_name + '.csv'
        df_data_process.to_csv(path_save_data_div, index=False)


def linear_regression(df_data, list_independent, str_dependent, str_descriptor_name='', weights=False,
                      display_plot=False):
    """ Linear regression.

    :param df_data: Pandas dataframe of the data.
    :param list_independent: List of the key names for the independent data in dataframe.
    :param str_dependent: String of the key name for the dependent data in dataframe.
    :param str_descriptor_name: String name added to statistics output.
    :param weights: True uses weights in dataframe and applies a weighted linear regression.
        False will pass a dataframe of the independent, dependent, and weights.
    :param display_plot: True will display the plot.
    :return Dataframe of the statistical results.
    """
    df_data_copy = df_data.copy()
    list_independent_copy = copy.deepcopy(list_independent)
    str_dependent_copy = str_dependent
    str_descriptor_name_copy = str_descriptor_name

    def isfloat(value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    df_empty_series = pd.Series([''])
    if weights:
        df_return_regression = pd.DataFrame({0: 'Linear Regression - WLS', 1: str_descriptor_name_copy}, index=[0])
    else:
        df_return_regression = pd.DataFrame({0: 'Linear Regression - OLS', 1: str_descriptor_name_copy}, index=[0])
    df_return_regression = pd.concat([df_return_regression, df_empty_series], ignore_index=True, sort=False)

    dict_data_independent = {}
    dict_data_dependent = {}

    for column in list_independent_copy:
        dict_data_independent[column] = copy.deepcopy(list(df_data_copy[column].dropna()))
    dict_data_dependent[str_dependent_copy] = copy.deepcopy(list(df_data_copy[str_dependent_copy].dropna()))

    list_weights = []
    if weights:
        list_weights = copy.deepcopy(list(df_data_copy['weights']))

    # Linear Regression
    df_x = pd.DataFrame(dict_data_independent, columns=list_independent_copy)
    df_y = pd.DataFrame(dict_data_dependent)

    df_x.insert(loc=0, column='const', value=1.0)

    if weights:
        sm_results = sm.WLS(df_y, df_x, weights=list_weights).fit()
    else:
        sm_results = sm.OLS(df_y, df_x).fit()

    sm_summary = sm_results.summary2()
    df_summary_0 = sm_summary.tables[0]
    df_return_regression = pd.concat([df_return_regression, df_summary_0, df_empty_series], ignore_index=True,
                                     sort=False)

    df_summary_1 = sm_summary.tables[1]
    df_summary_1 = df_summary_1.reset_index().T.reset_index().T.reset_index(drop=True)
    df_summary_1.at[0, 0] = ''

    df_return_regression = pd.concat([df_return_regression, df_summary_1, df_empty_series, df_empty_series],
                                     ignore_index=True,
                                     sort=False)

    sm_residuals = sm_results.resid
    sm_exogenous = sm_results.model.exog

    # Normality of the residuals
    # Shapiro-Wilk test
    df_shapiro_name = pd.DataFrame({0: 'Shapiro-Wilk test', 1: '(Normality of the residuals)'}, index=[0])
    df_sw_descriptor = pd.DataFrame({0: 'Null hypothesis:', 1: 'has normal distribution'}, index=[0])
    df_sw_labels = pd.DataFrame({0: 'W', 1: 'p-value'}, index=[0])
    sp_shapiro_test = stats.shapiro(sm_residuals)
    df_sw_results = pd.DataFrame(dict(lzip([0, 1], sp_shapiro_test)), index=[0])
    df_return_regression =\
        pd.concat([df_return_regression, df_shapiro_name, df_sw_descriptor, df_sw_labels, df_sw_results,
                   df_empty_series], ignore_index=True, sort=False)

    # Mutlicollinearity
    # Condition number
    df_conditionnum_name = pd.DataFrame({0: 'Condition number', 1: '(Mutlicollinearity)'}, index=[0])
    df_cn_descriptor = pd.DataFrame({0: 'If above 30:', 1: 'may have severe multicollinearity'}, index=[0])
    df_cn_labels = pd.DataFrame({0: 'number'}, index=[0])
    cn_results = np.linalg.cond(sm_exogenous[:, 1:])
    df_cn_results = pd.DataFrame({0: cn_results}, index=[0])
    df_return_regression = pd.concat([df_return_regression, df_conditionnum_name, df_cn_descriptor, df_cn_labels,
                                      df_cn_results, df_empty_series], ignore_index=True, sort=False)

    # Heteroskedasticity (constant variance)
    # White test
    df_white_name = pd.DataFrame({0: 'White test', 1: '(Heteroskedasticity)'}, index=[0])
    df_w_descriptor = pd.DataFrame({0: 'Null hypothesis:', 1: 'homoskedasticity'}, index=[0])
    df_w_labels = pd.DataFrame({0: 'Lagrange Multiplier statistic', 1: 'LM p-value', 2: 'F statistic', 3: 'F p-value'},
                               index=[0])
    w_results = sms.het_white(sm_residuals, sm_exogenous)
    df_w_results = pd.DataFrame(dict(lzip([0, 1, 2, 3], w_results)), index=[0])
    df_return_regression =\
        pd.concat([df_return_regression, df_white_name, df_w_descriptor, df_w_labels, df_w_results, df_empty_series],
                  ignore_index=True, sort=False)

    # Autocorrelation (independent error terms)
    # Durbin–Watson statistic
    df_durbinwatson_name = pd.DataFrame({0: 'Durbin-Watson test', 1: '(Autocorrelation)'}, index=[0])
    df_dw_descriptor = pd.DataFrame({0: 'Between 0-4', 1: '2 indicates no autocorrelation'}, index=[0])
    df_dw_labels = pd.DataFrame({0: 'test statistic'}, index=[0])
    dw_results = sms.durbin_watson(sm_residuals)
    df_dw_results = pd.DataFrame({0: dw_results}, index=[0])
    df_return_regression = pd.concat([df_return_regression, df_durbinwatson_name, df_dw_descriptor, df_dw_labels,
                                      df_dw_results, df_empty_series], ignore_index=True, sort=False)

    # Linearity (correct functional form)
    # Rainbow test
    df_rainbow_name = pd.DataFrame({0: 'Rainbow Test', 1: '(Linearity)'}, index=[0])
    df_rb_descriptor = pd.DataFrame({0: 'Null hypothesis:', 1: 'linear'}, index=[0])
    df_rb_labels = pd.DataFrame({0: 'F statistic', 1: 'p value'}, index=[0])
    rb_results = sms.linear_rainbow(sm_results)
    df_rb_results = pd.DataFrame(dict(lzip(df_rb_labels, rb_results)), index=[0])
    df_return_regression = pd.concat([df_return_regression, df_rainbow_name, df_rb_descriptor, df_rb_labels,
                                      df_rb_results, df_empty_series], ignore_index=True, sort=False)

    for col in df_return_regression:
        for num in range(len(df_return_regression[col])):
            item = df_return_regression[col][num]
            if isfloat(item):
                df_return_regression[col][num] = float(df_return_regression[col][num])

    if display_plot:
        list_keys_data_independent = list(dict_data_independent.keys())

        if len(list_keys_data_independent) > 1:
            str_xlabel = 'independent data'
        else:
            str_xlabel = list_keys_data_independent[0]
        str_ylabel = list(dict_data_dependent.keys())[0]
        plt_figure, plt_axis = plt.subplots()
        sm.graphics.plot_fit(sm_results, 1, ax=plt_axis, vlines=False)
        plt_axis.set_ylabel(str_ylabel)
        plt_axis.set_xlabel(str_xlabel)
        if str_descriptor_name_copy == '':
            plt_axis.set_title('Linear Regression')
        else:
            plt_axis.set_title('Linear Regression ' + str_descriptor_name_copy)
        plt.show()

    return df_return_regression


apply_regression(path_fasta=PATH_FASTA,
                 directory_load=DIRECTORY_LOAD,
                 file_load=FILE_LOAD,
                 orthodb_seq=ORTHODB_SEQ,
                 list_options=LIST_OPTIONS,
                 display_plot=DISPLAY_PLOT)
