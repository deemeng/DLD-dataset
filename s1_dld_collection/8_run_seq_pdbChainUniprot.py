import os
import pandas as pd

from utils.common import load_tab, dump_dict2json, read_json2dict
from loop.extract_loop import get_loops
from loop.dssp import generate_pdp_chain_pairs
from params import *
'''
!!!! note: loops including N-&C-terminus
Not useful currently

Apr 12, 2023
Di
'''
###
# 2. generate pdbid_chainid list
###
# {'1A04': ['A'], '1A0I': ['A'],...}
dict_pdb_chain_unp = read_json2dict(path_pdb_chain_unp)

# # ['1a04_a_uniprotID', '1a0i_a_uniprotID', '1a0p_a_uniprotID', ...]
list_pidChainUnp = []

for pid, dict_chain_unp in dict_pdb_chain_unp.items():
    pid = pid.lower()
    for chain, unps in dict_chain_unp.items():
        list_pidChainUnp += [f'{pid}_{chain.lower()}_{unp}' for unp in unps]
        
