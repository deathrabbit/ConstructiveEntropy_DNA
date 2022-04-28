# Process1_Retrieve.py
# Retrieves / updates NCBI taxonomy database
# and obtains Timetree.org phylogenetic tree files

__author__ = "Scott Matsumura"
__last_update__ = "4.15.2022"

import ssl
import os
import urllib3
from ete3 import NCBITaxa

PATH_FASTA = '/PyCharm/DNAEntropy/fasta/'
GET_NCBI = False
GET_TIMETREE = False


def retrieve_files(path_fasta, get_ncbi, get_timetree):
    """ Downloads NCBI taxonomy database and Timetree.org phylogenetic tree.

    :param path_fasta: Path to fasta directory. This will also be used to store files.
    :param get_ncbi: True will download / update NCBI taxonomy database
    :param get_timetree: True will download Timetree.org phylogenetic tree in Newick and NEXUS format
    """

    path_home = os.getenv('HOME')
    dir_timetree = 'zz_TimeTree'
    path_save = path_home + path_fasta + '/' + dir_timetree

    # Timetree of Life (2015) Newick format
    url_timetree2015_nwk = "http://www.timetree.org/public/data/TimetreeOfLife2015.nwk"
    save_nwk = 'TimetreeOfLife2015.nwk'

    # Timetree of Life (2015) NEXUS format
    url_timetree2015_tre = "http://www.timetree.org/public/data/TimetreeOfLife2015.tre"
    save_tre = 'TimetreeOfLife2015.tre'

    if not os.path.exists(path_save):
        os.makedirs(path_save)

    if get_ncbi:
        ssl_ctx = ssl.create_default_context()
        ssl_ctx.check_hostname = False
        ssl_ctx.verify_mode = ssl.CERT_NONE
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()
    if get_timetree:
        http = urllib3.PoolManager()

        response = http.request('GET', url_timetree2015_nwk)
        data = response.data.decode('utf-8')
        with open(path_save + save_nwk, 'wt') as f:
            f.write(data)
        f.close()

        response = http.request('GET', url_timetree2015_tre)
        data = response.data.decode('utf-8')
        with open(path_save + save_tre, 'wt') as f:
            f.write(data)
        f.close()


retrieve_files(path_fasta=PATH_FASTA, get_ncbi=GET_NCBI, get_timetree=GET_TIMETREE)
