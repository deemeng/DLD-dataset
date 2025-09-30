import os
import pandas as pd

from utils.common import load_tab, dump_dict2json, read_json2dict
from loop.extract_loop import get_loops
from loop.dssp import generate_pdp_chain_pairs
from params import *
'''
Aim to extract all loops (more than 4 continuous coil residues) from domain_linker/recollect_dsAll/dssp/full_chainid.
!!!! note: loops including N-&C-terminus

1st. smoothing SS elements with threshold=2 (<=2)
2nd. smoothing SS elements with threshold=3 (<=3)

Apr 12, 2023
Di
'''

###
# 1. Multi-domain file, from SCOP
###
# multi_domain = load_tab(path_tab_scop_FA_continuous_uni_multi_parse_domian)

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
    
###
# 3. find loops and save to a JSON file
###

print('extracting loops_smoothSS_first ...')
loops_smoothSS_first = get_loops(list_pidChainUnp, path_dssp_full_uniprot_chainid, ss_1=True, ss_2=False)
print('extracting loops_smoothSS_second ...')
loops_smoothSS_second = get_loops(list_pidChainUnp, path_dssp_full_uniprot_chainid, ss_1=False, ss_2=True)

print('saving ...')
dump_dict2json(loops_smoothSS_first, path_loops_all_smoothSS_first)
dump_dict2json(loops_smoothSS_second, path_loops_all_smoothSS_second)
print('Done')