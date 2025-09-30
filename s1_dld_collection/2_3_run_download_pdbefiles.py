import os
import numpy as np
import pandas as pd

import json
import wget

from params import *

from utils.common import load_tab
from loop.dssp import generate_pdp_chain_pairs
from loop.pdb_files import download_pdbe_updatedCif

'''
25 Apr, 2023
Di

Aim to download updated_mmcif files from pdbe.
e.g. 
    https://www.ebi.ac.uk/pdbe/entry-files/download/1a04_updated.cif
'''

multi_domain = load_tab(path_tab_scop_FA_continuous_uni_multi_parse_domian)
dict_pdb_chain = generate_pdp_chain_pairs(multi_domain)

# list of pdb entries
list_pdbid = list(dict_pdb_chain.keys())

download_pdbe_updatedCif(list_pdbid, path_pdbe_cif)