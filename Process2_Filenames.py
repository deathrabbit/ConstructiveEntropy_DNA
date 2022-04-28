# Process_2_Filenames.py
# Creates files that list the OrthoDB files to be processed and if correction is required
# based on a check with the NCBI database.

__author__ = "Scott Matsumura"
__last_update__ = "4.15.2022"

import os
import pickle
import bz2
import bisect
import re
from ete3 import NCBITaxa

PATH_FASTA = '/PyCharm/DNAEntropy/fasta/'

''' Select kingdom to process '''
KINGDOM = 'bacteria'
# KINGDOM = 'fungi'
# KINGDOM = 'metazoa'
# KINGDOM = 'plants'
# KINGDOM = 'protozoa'
# KINGDOM = 'viridae'


def process_filenames(kingdom, path_fasta):
    """ Indexes FASTA filenames
    Files are saved to "FileIndex" folder:

    "file_index_rawdata.txt": List of rawdata 'fs' file names saved to text file is same order as 'file_index_correct'.
    "file_index_correct.txt": List of NCBI corrected id saved to text file is same order as 'file_index_rawdata
    "file_index_correct_taxid.txt": List of NCBI corrected id without alternates
    "file_id_corrections.txt": List of all NCBI corrected changes saved to text file
    "file_correct_dict.pbz2":
        key: NCBI taxid
        value: List of corrected file names for alternative sequencing
    "file_counts_alternates.txt": List of counts of alternative sequencing
    "file_not_found.txt": List of OrthoDB species files not found in NCBI database.
        These files will not be included and further processed.

    :param kingdom: kingdom name to be processed.
    :param path_fasta: location of fasta directory from the home directory, which will also be used as save directory.
    """

    ''' Nested functions '''
    def save_bz2(filename, save_file):
        """ Saves file to binary file with bzip2 compression

        :param filename: file name
        :param save_file: file to be saved
        """
        with bz2.BZ2File(filename=filename, mode='wb', compresslevel=9) as f:
            pickle.dump(save_file, f)
        f.close()

    def save_list_to_text(filename, save_list):
        """ Saves array list to text file

        :param filename: file name
        :param save_list: array list to be saved
        """
        with open(file=filename, mode='wt', encoding='utf-8') as f:
            f.write('\n'.join(save_list))
        f.close()

    ''' Start of "def process_filenames '''

    # Directory and file names
    dir_rawdata = 'RawData'
    dir_fileindex = 'FileIndex'
    ext_rawdata = 'fs'
    ext_compression = 'bz2'

    ncbi = NCBITaxa()

    path_home = os.getenv('HOME')
    path_load = path_home + path_fasta + kingdom + '/' + dir_rawdata + '/'
    path_index = path_home + path_fasta + kingdom + '/' + dir_fileindex + '/'

    ext_rawdata = ext_rawdata
    ext_compression = ext_compression

    file_index_rawdata_txt = 'file_index_rawdata.txt'
    file_index_correct_txt = 'file_index_correct.txt'
    file_index_correct_taxid_txt = 'file_index_correct_taxid.txt'
    file_id_corrections_txt = 'file_id_corrections.txt'
    file_correct_dict_bz2 = 'file_correct_dict.bz2'
    file_counts_alternates_txt = 'file_counts_alternates.txt'
    file_not_found_txt = 'file_not_found.txt'

    print("Processing file names.")

    list_index_rawdata = sorted([f for f in os.listdir(path_load) if f.endswith(ext_rawdata)])

    if not os.path.exists(path_index):
        os.makedirs(path_index)

    list_index_correct = []
    list_index_taxid = []
    list_species_name_old = []
    list_species_name_new = []
    dict_index_correct = {}
    dict_alternates_count = {}
    list_not_found = []
    delimiters = '.', '_'
    regex_pattern = '|'.join(map(re.escape, delimiters))
    for index_file_old in list_index_rawdata:
        line_list = re.split(regex_pattern, index_file_old)
        index_key = list(ncbi.get_taxid_translator([line_list[0]]).keys())
        if len(index_key) == 0:
            list_not_found.append(index_file_old)
        else:
            taxid_name = ncbi.get_taxid_translator([line_list[0]])[int(line_list[0])]
            taxid_num = str((ncbi.get_name_translator([taxid_name])[taxid_name])[0])
            index_file_new = index_file_old
            if len(line_list) == 2:
                index_file_new = taxid_num + '_0.' + ext_compression
                if dict_alternates_count.get('0') is None:
                    dict_alternates_count['0'] = 1
                else:
                    dict_alternates_count['0'] = dict_alternates_count['0'] + 1
            if len(line_list) == 3:
                index_file_new = taxid_num + '_' + line_list[1] + '.' + ext_compression
                if dict_alternates_count.get(line_list[1]) is None:
                    dict_alternates_count[line_list[1]] = 1
                else:
                    dict_alternates_count[line_list[1]] = dict_alternates_count[line_list[1]] + 1

            if line_list[0] != taxid_num:
                list_species_name_old.append(index_file_old)
                list_species_name_new.append(index_file_new)
            list_index_correct.append(index_file_new)

            if dict_index_correct.get(taxid_num) is None:
                dict_index_correct[taxid_num] = [index_file_new]
                list_index_taxid.append(taxid_num)
            else:
                bisect.insort(dict_index_correct[taxid_num], index_file_new)

    for index_not_found in list_not_found:
        list_index_rawdata.remove(index_not_found)

    save_path_list_index_rawdata = path_index + file_index_rawdata_txt
    save_list_to_text(save_path_list_index_rawdata, list_index_rawdata)

    save_path_list_index_correct = path_index + file_index_correct_txt
    save_list_to_text(save_path_list_index_correct, list_index_correct)

    save_path_dict_index_correct = path_index + file_correct_dict_bz2
    save_bz2(save_path_dict_index_correct, dict_index_correct)

    save_path_list_index_taxid = path_index + file_index_correct_taxid_txt
    save_list_to_text(save_path_list_index_taxid, list_index_taxid)

    list_id_changed_to = ['Corrected NCBI ID', 'Old -> New', '']
    for index_old, index_new in zip(list_species_name_old, list_species_name_new):
        list_id_changed_to.append(index_old + ' -> ' + index_new)
    save_path_list_id_changed_to = path_index + file_id_corrections_txt
    save_list_to_text(save_path_list_id_changed_to, list_id_changed_to)

    counts_alternates_items = list(dict_alternates_count.items())
    list_counts_alternates = ['Number of multiple species by name extension', '']
    for counts in counts_alternates_items:
        list_counts_alternates.append('_' + str(counts[0] + ' : ' + str(counts[1])))
    save_path_list_counts_alternates = path_index + file_counts_alternates_txt
    save_list_to_text(save_path_list_counts_alternates, list_counts_alternates)

    list_not_found_species = ['List of OrthoDB species files not found in NCBI database', '']
    for str_not_fount in list_not_found:
        list_not_found_species.append(str_not_fount)
    save_path_list_not_found_species = path_index + file_not_found_txt
    save_list_to_text(save_path_list_not_found_species, list_not_found_species)

    print("Done.")


process_filenames(kingdom=KINGDOM, path_fasta=PATH_FASTA)
