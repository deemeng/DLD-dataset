import os
import pandas as pd

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBList

from loop.pdb_files import get_PDBfiles
from loop.dssp import save_SS
from utils.common import load_tab, save_tab
from params import *

'''
14 Mar, 2023
Di

Aim to save DSSP information of each PDB entry to an JSON file.

1. Secondary structure
2. Hbonds information
3. Beta-bridge information

'''

###
# Save all secondary structures to files. name as pdbid.json
# 1. generate dicts for each PDBID
# 2. save them to json files
###

list_pdbid = load_tab(path_tab_scop_FA_continuous_uni_multi_domian_lenResolution_filtered)['FA-PDBID'].unique() # path_tab_scop_FA_continuous_uni_multi_domian
error_list = save_SS(list_pdbid[:2], path_pdb_cif, path_dssp_pdbid)
# save the list of dicts to json files
# dump_dicts2jsons(list_ss, folder=path_dssp_pdbid, chain=False)

# some pdbID cannot be correctly parsed by bio.PDB.DSSP
with open(path_dssp_error, 'w') as f:
    for line in error_list:
        f.write(f"{line[0]}, {line[1]}\n")
