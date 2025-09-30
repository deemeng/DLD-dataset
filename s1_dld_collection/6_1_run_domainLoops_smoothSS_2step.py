import os
import pandas as pd
import numpy as np
import json

from loop.dssp import save_chains, generate_pdp_chain_pairs
from utils.common import load_tab, dump_dict2json, int_domainREG, save_tab, read_json2dict
from loop.extract_loop import get_loops
from domain_linker.extract_linkers import get_domain_loops, save_missingPercentage

from params import *

'''
Aim to 
1. extract & save inter-domain loops
2. calculate & save inter-domain loops with missing percentage
'''

multi_domain = load_tab(path_tab_scop_FA_continuous_uni_multi_parse_domian)
with open(path_loops_all_smoothSS_second, 'r') as f:
    loops = json.load(f)
    
# domain_loops, adjacent_domains, num_loops 
domain_loops, _, remove_pdbchain, num_loops = get_domain_loops(multi_domain, loops)
print(f'The number of loops: {num_loops}')
print('Saving inter-domain loops ...')
dump_dict2json(domain_loops, path_loops_domain_smoothSS_second)

# calculate & save inter-domain loops with missing percentage
print('calculate & save inter-domain loops with missing percentage ...')
save_missingPercentage(domain_loops, path_dssp_full_uniprot_chainid, path_loops_domain_smoothSS_second_perMiss)
print('Done!')

print(f'SIFT & SCOP does not match: {remove_pdbchain}')
print('Update path_tab_scop_FA_continuous_uni_multi_parse_domian')
# ['6S9U']
for pdbchain in remove_pdbchain:
    pdb = pdbchain[0]
    chain = pdbchain[1]
    multi_domain = multi_domain[~((multi_domain['FA-PDBID']==pdb)&(multi_domain['FA-CHAINID']==chain))]
# save
save_tab(multi_domain, path_tab_scop_FA_continuous_uni_multi_parse_domian)

# Update path_pdb_chain_unp
print('Update path_pdb_chain_unp')
dict_pdb_chain_unp = {}
df_dict = multi_domain[['FA-PDBID', 'FA-CHAINID', 'FA-UNIID']]
for pdbid, sub_df in df_dict.groupby("FA-PDBID"):
    pdbid_lower = pdbid.lower()
    chain_dict = {}
    for chain_id, chain_df in sub_df.groupby("FA-CHAINID"):
        unique_uniprot = sorted(set(chain_df["FA-UNIID"]))  # remove duplicates
        chain_dict[chain_id] = unique_uniprot
    dict_pdb_chain_unp[pdbid_lower] = chain_dict

# save dict_pdb_chain_unp
dump_dict2json(dict_pdb_chain_unp, path_pdb_chain_unp)
print('Done!')