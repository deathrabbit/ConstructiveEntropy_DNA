# Process_4_Check.py
# Produces three CSV files to help check for errors:
#     (kingdom)_amino.csv is a count of amino acids.
#     (kingdom)_unique.csv is a count of the repeated proteins by number of repeats.
#     (kingdom)_entropy.csv is the amino, protein, and constructive entropy values plus divergence time if available.

__author__ = "Scott Matsumura"
__last_update__ = "4.15.2022"

import os
import pickle
import bz2
import pandas as pd

PATH_FASTA = '/My Drive/PyCharm/DNAEntropy/fasta/'

""" Select kingdom to process"""
KINGDOM = 'bacteria'
# KINGDOM = 'fungi'
# KINGDOM = 'metazoa'
# KINGDOM = 'plants'
# KINGDOM = 'protozoa'
# KINGDOM = 'viridae'


def csv_check(kingdom, path_fasta):
    """ Convert data into csv graphs
    Creates CSV files to check the correctness of the processed species files.
    (kingdom)_amino.csv is a count of amino acids.
    (kingdom)_unique.csv is a count of the repeated proteins by number of repeats.
    (kingdom)_entropy.csv is the amino, protein, and constructive entropy values plus divergence time if available.

    :param kingdom: kingdom name to be processed.
    :param path_fasta: location of fasta directory from the home directory, which will also be used as save directory.
    """

    ''' Nested functions '''
    def load_bz2(filename):
        """ Loads binary bzip2 file to array list

        :param filename: file name
        :return: file
        """
        with bz2.BZ2File(filename=filename, mode='rb') as f:
            return_file = pickle.load(f)
        f.close()
        return return_file

    def load_text_to_list(filename):
        """ Loads text file to array list

        :param filename: file name
        :return: array list
        """
        with open(file=filename, mode='rt') as f:
            return_array = f.read().splitlines()
        f.close()
        return return_array

    ''' Start of "def process_species '''

    # Directory and file names
    dir_process_species = 'ProcessSpecies'
    dir_fileindex = 'FileIndex'
    dir_csvcheck = 'zz_CSVcheck'

    path_home = os.getenv('HOME')
    path_index = path_home + path_fasta + kingdom + '/' + dir_fileindex + '/'
    path_save = path_home + path_fasta + '/' + dir_csvcheck + '/'
    path_load = path_home + path_fasta + kingdom + '/' + dir_process_species + '/'

    file_amino_csv = kingdom + '_amino_count.csv'
    file_unique_csv = kingdom + '_unique_count.csv'
    file_entropy_csv = kingdom + '_entropy_count.csv'
    file_index_correct_txt = 'file_index_correct.txt'
    list_index_corrected = load_text_to_list(path_index + file_index_correct_txt)

    if not os.path.exists(path_save):
        os.makedirs(path_save)

    print("Processing CSV Check:")
    print("Loading " + kingdom + " files...")

    str_column_name = ['Species', 'OrthoDB_name']
    str_column_amino = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    str_column_entropy = ['amino_entropy', 'protein_entropy', 'constructive_entropy', 'divergence_time']
    df_name = pd.DataFrame(columns=str_column_name)
    df_amino = pd.DataFrame(columns=str_column_amino)
    df_unique = pd.DataFrame()
    df_entropy = pd.DataFrame(columns=str_column_entropy)
    size_index_corrected = len(list_index_corrected)

    # Retrieve data from each file
    for int_index, str_file in enumerate(list_index_corrected, 1):
        path_load_species_bz2 = path_load + str_file
        dict_species_data = load_bz2(path_load_species_bz2)
        dict_data = dict_species_data['data']
        dict_line_amino = dict_species_data['amino']
        dict_line_name = {'Species': [dict_data['species_name']], 'OrthoDB_name': [dict_data['orthodb_id']]}
        dict_line_entropy = {'amino_entropy': [dict_data['amino_entropy']],
                             'protein_entropy': [dict_data['protein_entropy']],
                             'constructive_entropy': [dict_data['constructive_entropy']],
                             'divergence_time': [dict_data['divergence_time']]}
        dict_line_unique = dict_species_data['unique']

        # Wraps all values into lists to concatenate into a DataFrame
        for str_key in dict_line_amino.keys():
            dict_line_amino[str_key] = [dict_line_amino[str_key]]
        for str_key in dict_line_unique.keys():
            dict_line_unique[str_key] = [dict_line_unique[str_key]]

        # Concatenate dictionary into DataFrame
        # Process is slow because the concatenation deep copies each time.
        # It is done this way to align the values correctly.
        df_name = pd.concat([df_name, pd.DataFrame.from_dict(dict_line_name)], ignore_index=True)
        df_amino = pd.concat([df_amino, pd.DataFrame.from_dict(dict_line_amino)], ignore_index=True)
        df_unique = pd.concat([df_unique, pd.DataFrame.from_dict(dict_line_unique)], ignore_index=True, sort=True)
        df_entropy = pd.concat([df_entropy, pd.DataFrame.from_dict(dict_line_entropy)], ignore_index=True)
        print(str_file + ' : ' + str(int_index) + ' of ' + str(size_index_corrected))

    # Save to file
    df_amino[list(df_amino)] = df_amino[list(df_amino)].fillna(0.0).astype(int)
    df_amino = pd.concat([df_name, df_amino], axis=1)
    path_save_amino = path_save + file_amino_csv
    df_amino.to_csv(path_save_amino, index=False)

    df_unique[list(df_unique)] = df_unique[list(df_unique)].fillna(0.0).astype(int)
    df_unique = pd.concat([df_name, df_unique], axis=1)
    path_save_unique = path_save + file_unique_csv
    df_unique.to_csv(path_save_unique, index=False)

    df_entropy = pd.concat([df_name, df_entropy], axis=1)
    path_save_entropy = path_save + file_entropy_csv
    df_entropy.to_csv(path_save_entropy, index=False)

    print("Done.")


csv_check(kingdom=KINGDOM, path_fasta=PATH_FASTA)
