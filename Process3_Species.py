# Process_3_Species.py
# To be used with OrthoDB v10
# Processes each species OrthoDB files and saves the data

__author__ = "Scott Matsumura"
__last_update__ = "4.15.2022"

import os
import re
import pickle
import bz2
from collections import Counter
from Bio import SeqIO, Phylo
from ete3 import NCBITaxa

PATH_FASTA = '/My Drive/PyCharm/DNAEntropy/fasta/'

""" Select kingdom to process"""
KINGDOM = 'bacteria'
CLEAN_CHARACTERS = False
CLEAN_LOWERCASE = False

# KINGDOM = 'fungi'
# CLEAN_CHARACTERS = False
# CLEAN_LOWERCASE = False

# KINGDOM = 'metazoa'
# CLEAN_CHARACTERS = True
# CLEAN_LOWERCASE = False

# KINGDOM = 'plants'
# CLEAN_CHARACTERS = True
# CLEAN_LOWERCASE = False

# KINGDOM = 'protozoa'
# CLEAN_CHARACTERS = True
# CLEAN_LOWERCASE = False

# KINGDOM = 'viridae'
# CLEAN_CHARACTERS = False
# CLEAN_LOWERCASE = False


def process_species(kingdom, path_fasta, clean_characters=False, clean_lowercase=False):
    """ Loads FASTA files then processes and saves them into Python bzip2 files
    File is saved as a dictionary with the following keys.
    title:
        key: protein
        value: list of titles of that particular protein
    protein:
        key: protein
        value: protein count
    amino:
        key: amino acid letter
        value: amino acid count
    unique:
        key: count of repeated proteins
        value: number of repeated proteins for that count
    data:
        key: "amino_count", value: total number of amino acids including duplicates
        key: "amino_count_unique", value: number of amino acids excluding duplicates
        key: "amino_edge_count", value: number of amino differences (graph edges)
        key: "protein_count", value: total number of proteins including duplicates
        key: "protein_count_unique", value: total number of excluding duplicates
        key: "protein_edge_count", value: number of protein differences (graph edges)
        key: "amino_entropy", value: amino acid entropy
        key: "protein_entropy", value: protein entropy
        key: "constructive_entropy", value: constructive entropy
        key: "divergence_time", value: divergence time if available, 'None' otherwise
        key: "orthodb_id", value: string of OrthoDB species name
        key: "sequence_num", value: int of OrthoDB sequence number
        key: "species_id", value: int of species NCBI id
        key: "species_name", value: string of species name
        key: "lineage_id", value: list of lineage by NCBI id
        key: "lineage_name", value: list of lineage by name

    :param kingdom: kingdom name to be processed.
    :param path_fasta: location of fasta directory from the home directory, which will also be used as save directory.
    :param clean_characters: Removes non alphas from protein. Default is false.
    :param clean_lowercase: Replaces lowercase amino acid by 'X'. Default is false.
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

    def load_text_to_list(filename):
        """ Loads text file to array list
        :param filename: file name
        :return: array list
        """
        with open(file=filename, mode='rt') as f:
            return_array = f.read().splitlines()
        f.close()
        return return_array

    def clean_characters_protein(line):
        error_symbol = set()
        for character in line:
            if not character.isalpha():
                error_symbol.add(character)
        for error in error_symbol:
            line = line.replace(error, '')
        return line

    def clean_lowercase_protein(line):
        char_array = list(line)
        for x in range(len(char_array)):
            if char_array[x].islower():
                char_array[x] = 'X'
        line = ''.join(char_array)
        return line

    def difference_entropy_by_counts(list_proportions):
        """ Calculates difference entropy by counts
        :param list_proportions: List of element counts
        :return: int_total (N), int_unique (M), int_edge (|E(G)|), flt_entropy (H)
        """
        int_total = sum(list_proportions)
        int_unique = len(list_proportions)
        int_total_2 = int_total * int_total
        int_edge = (int_total_2 / 2) - (sum(x*x for x in list_proportions) / 2)
        flt_entropy = (2 * int_edge) / int_total_2
        return int_total, int_unique, int_edge, flt_entropy

    ''' Start of "def process_species '''

    # Directory and file names
    dir_rawdata = 'RawData'
    dir_fileindex = 'FileIndex'
    dir_process_species = 'ProcessSpecies'
    dir_timetree = 'zz_TimeTree'

    path_home = os.getenv('HOME')
    path_load = path_home + path_fasta + kingdom + '/' + dir_rawdata + '/'
    path_save = path_home + path_fasta + kingdom + '/' + dir_process_species + '/'
    path_index = path_home + path_fasta + kingdom + '/' + dir_fileindex + '/'
    path_timetree_nwk = path_home + path_fasta + '/' + dir_timetree + '/'

    file_index_rawdata_txt = 'file_index_rawdata.txt'
    file_index_correct_txt = 'file_index_correct.txt'
    file_timetree_nwk = 'TimetreeOfLife2015.nwk'

    list_index_rawdata = load_text_to_list(path_index + file_index_rawdata_txt)
    list_index_corrected = load_text_to_list(path_index + file_index_correct_txt)
    tree_nwk = Phylo.read(path_timetree_nwk + file_timetree_nwk, 'newick')

    int_index_size = len(list_index_rawdata)
    int_divergence_count = 0

    ncbi = NCBITaxa()

    if not os.path.exists(path_save):
        os.makedirs(path_save)

    print("Processing species:")
    for int_index, (str_index_file_old, str_index_file_new)\
            in enumerate(zip(list_index_rawdata, list_index_corrected), 1):
        str_orthodb_id, str_extension_new = str_index_file_new.split('.')
        dict_title = {}
        counter_protein = Counter()
        counter_amino = Counter()
        counter_unique = Counter()
        with open(file=path_load + str_index_file_old, mode='rt') as file:
            for record in SeqIO.parse(file, 'fasta'):
                str_seq = str(record.seq)
                str_id_title = str(record.id)
                id_title_name, id_title_index = str_id_title.split(':')
                str_id_title = str_orthodb_id + ':' + id_title_index
                if clean_characters:
                    str_seq = clean_characters_protein(str_seq)
                if clean_lowercase:
                    str_seq = clean_lowercase_protein(str_seq)
                counter_protein.update([str_seq])
                counter_amino.update(str_seq)
                if record.id in dict_title:
                    dict_title[str_seq].append(str_id_title)
                else:
                    dict_title[str_seq] = [str_id_title]
        file.close()

        list_protein_values = list(counter_protein.values())
        int_protein_total, int_protein_count_unique, int_protein_edge, flt_protein_entropy = \
            difference_entropy_by_counts(list_protein_values)

        list_amino_values = list(counter_amino.values())
        int_amino_count_total, int_amino_count_unique, int_amino_edge, flt_amino_entropy = \
            difference_entropy_by_counts(list_amino_values)

        flt_constructive_entropy = flt_protein_entropy - flt_amino_entropy

        counter_unique.update(list_protein_values)

        dict_protein = dict(counter_protein)
        dict_amino = dict(counter_amino)
        dict_unique = dict(counter_unique)

        str_taxid, str_sequence_num = str_orthodb_id.split('_')
        int_taxid = int(str_taxid)
        int_sequence_num = int(str_sequence_num)
        str_species_name_ncbi = list(ncbi.get_taxid_translator([int_taxid]).values())[0]
        list_lineage_ncbi_id = ncbi.get_lineage(int_taxid)
        list_lineage_ncbi_name = []
        for n in range(1, len(list_lineage_ncbi_id)):
            lineage_names_ncbi = list(ncbi.get_taxid_translator([list_lineage_ncbi_id[n]]).values())[0]
            list_lineage_ncbi_name.append(lineage_names_ncbi)

        str_species_name_ncbi_nwk = re.sub(r'([^\s\w]|_)+', '', str_species_name_ncbi)
        str_species_name_ncbi_nwk = str_species_name_ncbi_nwk.replace(' ', '_')

        tree_clades = tree_nwk.find_clades(name=str_species_name_ncbi_nwk)
        flt_divergence_time = None
        for clade in tree_clades:
            flt_divergence_time = float(list(clade.depths().values())[0])
            int_divergence_count += 1

        data = {'amino_count': int_amino_count_total,
                'amino_count_unique': int_amino_count_unique,
                'amino_edge_count': int_amino_edge,
                'protein_count': int_protein_total,
                'protein_count_unique': int_protein_count_unique,
                'protein_edge_count': int_protein_edge,
                'amino_entropy': flt_amino_entropy,
                'protein_entropy': flt_protein_entropy,
                'constructive_entropy': flt_constructive_entropy,
                'divergence_time': flt_divergence_time,
                'orthodb_id': str_orthodb_id,
                'sequence_num': int_sequence_num,
                'species_id': int_taxid,
                'species_name': str_species_name_ncbi,
                'lineage_id': list_lineage_ncbi_id,
                'lineage_name': list_lineage_ncbi_name}

        save_path_dict_species_data = path_save + str_index_file_new
        dict_species_data =\
            {'title': dict_title, 'protein': dict_protein, 'amino': dict_amino, 'unique': dict_unique, 'data': data}
        save_bz2(save_path_dict_species_data, dict_species_data)
        print("Processing: " + str_index_file_old + " (" + "div_ct " + str(int_divergence_count) + ")" + " : " +
              str(int_index) + " of " + str(int_index_size))
    print("Done.")


process_species(kingdom=KINGDOM, path_fasta=PATH_FASTA, clean_characters=CLEAN_CHARACTERS,
                clean_lowercase=CLEAN_LOWERCASE)
