# Process_5_CompileData.py

__author__ = "Scott Matsumura"
__last_update__ = "4.15.2022"

import os
import pickle
import bz2
import pandas as pd
from ete3 import NCBITaxa

PATH_FASTA = '/My Drive/PyCharm/DNAEntropy/fasta/'
LIST_KINGDOM = ['bacteria', 'fungi', 'metazoa', 'plants', 'protozoa', 'viridae']


def compile_to_csv(list_kingdom, path_fasta):
    """ Convert data into a single csv graph to be used to run statistical analysis on.

    :param list_kingdom: List of kingdom names to be processed.
    :param path_fasta: Location of fasta directory from the home directory, which will also be used as save directory.
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
    dir_csvdata = 'zz_CSVdata'
    file_csvdata = 'entropy_data.csv'
    file_index_correct_txt = 'file_index_correct.txt'

    path_home = os.getenv('HOME')
    path_save = path_home + path_fasta + '/' + dir_csvdata + '/'

    if not os.path.exists(path_save):
        os.makedirs(path_save)

    ncbi = NCBITaxa()

    # Column names of CSV file
    list_columns = ['NCBI_id',                  # 0
                    'OrthoDB_seq',              # 1
                    'species_name',             # 2
                    'amino_entropy',            # 3
                    'protein_entropy',          # 4
                    'constructive_entropy',     # 5
                    'divergence_time',          # 6
                    'amino_edge_E(G)',          # 7
                    'amino_total_N',            # 8
                    'protein_edge_E(G)',        # 9
                    'protein_total_N',          # 10
                    'superkingdom',             # 11
                    'superkingdom_name',        # 12
                    'kingdom',                  # 13
                    'kingdom_name',             # 14
                    'subkingdom',               # 15
                    'subkingdom_name',          # 16
                    'phylum',                   # 17
                    'phylum_name',              # 18
                    'subphylum',                # 19
                    'subphylum_name',           # 20
                    'superclass',               # 21
                    'superclass_name',          # 22
                    'class',                    # 23
                    'class_name',               # 24
                    'subclass',                 # 25
                    'subclass_name',            # 26
                    'order',                    # 27
                    'order_name',               # 28
                    'family',                   # 29
                    'family_name',              # 30
                    'genus',                    # 31
                    'genus_name']               # 32
    df_data_compile = pd.DataFrame(columns=list_columns)

    print("Processing Compile to CSV:")
    for str_kingdom in list_kingdom:
        path_index = path_home + path_fasta + str_kingdom + '/' + dir_fileindex + '/'
        path_load = path_home + path_fasta + str_kingdom + '/' + dir_process_species + '/'
        list_index_corrected = load_text_to_list(path_index + file_index_correct_txt)

        size_index_corrected = len(list_index_corrected)
        for int_index, str_file in enumerate(list_index_corrected, 1):
            path_load_species_bz2 = path_load + str_file
            dict_species_data = load_bz2(path_load_species_bz2)
            dict_data = dict_species_data['data']
            list_lineage = dict_data['lineage_id']
            dict_rank = ncbi.get_rank(list_lineage)
            dict_rank = {value: key for key, value in dict_rank.items()}
            list_lineage_id = [dict_rank.get('superkingdom'),
                               dict_rank.get('kingdom'),
                               dict_rank.get('subkingdom'),
                               dict_rank.get('phylum'),
                               dict_rank.get('subphylum'),
                               dict_rank.get('superclass'),
                               dict_rank.get('class'),
                               dict_rank.get('subclass'),
                               dict_rank.get('order'),
                               dict_rank.get('family'),
                               dict_rank.get('genus')]
            list_lineage_id = [1 if x is None else x for x in list_lineage_id]
            dict_lineage_name = ncbi.get_taxid_translator(list_lineage_id)
            list_lineage_name = [dict_lineage_name[list_lineage_id[0]],
                                 dict_lineage_name[list_lineage_id[1]],
                                 dict_lineage_name[list_lineage_id[2]],
                                 dict_lineage_name[list_lineage_id[3]],
                                 dict_lineage_name[list_lineage_id[4]],
                                 dict_lineage_name[list_lineage_id[5]],
                                 dict_lineage_name[list_lineage_id[6]],
                                 dict_lineage_name[list_lineage_id[7]],
                                 dict_lineage_name[list_lineage_id[8]],
                                 dict_lineage_name[list_lineage_id[9]],
                                 dict_lineage_name[list_lineage_id[10]]]
            list_lineage_id = [x if x != 1 else None for x in list_lineage_id]
            list_lineage_name = [x if x != 'root' else None for x in list_lineage_name]

            # dict_line = {list_columns[0]: dict_data['species_id'],
            #              list_columns[1]: dict_data['sequence_num'],
            #              list_columns[2]: dict_data['species_name'],
            #              list_columns[3]: dict_data['amino_entropy'],
            #              list_columns[4]: dict_data['protein_entropy'],
            #              list_columns[5]: dict_data['constructive_entropy'],
            #              list_columns[6]: dict_data['divergence_time'],
            #              list_columns[7]: dict_data['amino_edge_count'],
            #              list_columns[8]: dict_data['amino_count'],
            #              list_columns[9]: dict_data['protein_edge_count'],
            #              list_columns[10]: dict_data['protein_count'],
            #              list_columns[11]: list_lineage_id[0],
            #              list_columns[12]: list_lineage_name[0],
            #              list_columns[13]: list_lineage_id[1],
            #              list_columns[14]: list_lineage_name[1],
            #              list_columns[15]: list_lineage_id[2],
            #              list_columns[16]: list_lineage_name[2],
            #              list_columns[17]: list_lineage_id[3],
            #              list_columns[18]: list_lineage_name[3],
            #              list_columns[19]: list_lineage_id[4],
            #              list_columns[20]: list_lineage_name[4],
            #              list_columns[21]: list_lineage_id[5],
            #              list_columns[22]: list_lineage_name[5],
            #              list_columns[23]: list_lineage_id[6],
            #              list_columns[24]: list_lineage_name[6],
            #              list_columns[25]: list_lineage_id[7],
            #              list_columns[26]: list_lineage_name[7],
            #              list_columns[27]: list_lineage_id[8],
            #              list_columns[28]: list_lineage_name[8],
            #              list_columns[29]: list_lineage_id[9],
            #              list_columns[30]: list_lineage_name[9],
            #              list_columns[31]: list_lineage_id[10],
            #              list_columns[32]: list_lineage_name[10]}
            # df_data_compile = df_data_compile.append(dict_line, ignore_index=True)

            dict_line = {list_columns[0]: [dict_data['species_id']],
                         list_columns[1]: [dict_data['sequence_num']],
                         list_columns[2]: [dict_data['species_name']],
                         list_columns[3]: [dict_data['amino_entropy']],
                         list_columns[4]: [dict_data['protein_entropy']],
                         list_columns[5]: [dict_data['constructive_entropy']],
                         list_columns[6]: [dict_data['divergence_time']],
                         list_columns[7]: [dict_data['amino_edge_count']],
                         list_columns[8]: [dict_data['amino_count']],
                         list_columns[9]: [dict_data['protein_edge_count']],
                         list_columns[10]: [dict_data['protein_count']],
                         list_columns[11]: [list_lineage_id[0]],
                         list_columns[12]: [list_lineage_name[0]],
                         list_columns[13]: [list_lineage_id[1]],
                         list_columns[14]: [list_lineage_name[1]],
                         list_columns[15]: [list_lineage_id[2]],
                         list_columns[16]: [list_lineage_name[2]],
                         list_columns[17]: [list_lineage_id[3]],
                         list_columns[18]: [list_lineage_name[3]],
                         list_columns[19]: [list_lineage_id[4]],
                         list_columns[20]: [list_lineage_name[4]],
                         list_columns[21]: [list_lineage_id[5]],
                         list_columns[22]: [list_lineage_name[5]],
                         list_columns[23]: [list_lineage_id[6]],
                         list_columns[24]: [list_lineage_name[6]],
                         list_columns[25]: [list_lineage_id[7]],
                         list_columns[26]: [list_lineage_name[7]],
                         list_columns[27]: [list_lineage_id[8]],
                         list_columns[28]: [list_lineage_name[8]],
                         list_columns[29]: [list_lineage_id[9]],
                         list_columns[30]: [list_lineage_name[9]],
                         list_columns[31]: [list_lineage_id[10]],
                         list_columns[32]: [list_lineage_name[10]]}
            df_data_compile = pd.concat([df_data_compile, pd.DataFrame.from_dict(dict_line)], ignore_index=True)

            print(str_file + ' : ' + str(int_index) + ' of ' + str(size_index_corrected) + ' ' + str_kingdom)

        path_save_compile = path_save + file_csvdata
        df_data_compile.to_csv(path_save_compile, index=False)

    print("Done.")


compile_to_csv(list_kingdom=LIST_KINGDOM, path_fasta=PATH_FASTA)
