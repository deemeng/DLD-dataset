import os
import pandas as pd

from loop.dssp import save_chains, generate_pdp_chain_pairs
from utils.common import load_tab, lower_pdbid
from params import *
'''
14 Mar, 2023
Di

!!!! similar to run_saveSSchianid.py, only the file paths are different.
Aim to split chains in a PDB entry. save each chain to a JSON file.
- Only chains wilth multiple domains are saved. 
- chains with one or zero domain are ignored.
'''

###
# 1. Load files & get pdbids
###
multi_domain = load_tab(path_tab_scop_FA_continuous_uni_multi_parse_domian)
list_PDBID = multi_domain['FA-PDBID'].unique()
list_pdbid = lower_pdbid(list_PDBID)

###
# 2. process and save JSON files
###
dict_pdb_chain = generate_pdp_chain_pairs(multi_domain)

print(f'Saving ... ')

for pdbid in list_pdbid:
    try:
        # Load json file to a DataFrame
        pdbid = pdbid.lower()
        ss_path = os.path.join(path_dssp_pdbid, pdbid + '.json')
        df_ss = pd.read_json(ss_path)
        save_chains(pdbid, df_ss, path_dssp_chainid, dict_pdb_chain)
    except Exception as e:
        print(f'{pdbid}: {e}')
