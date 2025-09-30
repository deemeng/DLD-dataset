import os
import pandas as pd

from loop.dssp import save_chains_full, generate_pdp_chain_pairs
from utils.common import load_tab, lower_pdbid, save_tab, dump_dict2json
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
print(f'Number of PDB entries: {len(list_PDBID)}')
list_pdbid = lower_pdbid(list_PDBID)

###
# 2. process and save JSON files
###
dict_pdb_chain = generate_pdp_chain_pairs(multi_domain)

print(f'Saving ... ')

list_pdbChainID_no_unp = []
dict_pdb_chain_unp = {}

for pdbid in list_pdbid:
    # try:
    # Load json file to a DataFrame
    pdbid = pdbid.lower()
    print(pdbid)
    ss_path = os.path.join(path_dssp_full_pdbid, pdbid + '.json')
    df_ss = pd.read_json(ss_path)
    
    # {'list_pdbChainID_no_unp': list_pdbChainID_no_unp, 'chain_unpID': chain_unpID}
    result = save_chains_full(pdbid, df_ss, path_dssp_full_uniprot_chainid, dict_pdb_chain, path_pdbe_uniprot_pdbid)
    list_pdbChainID_no_unp += result['list_pdbChainID_no_unp'] 
    dict_pdb_chain_unp[pdbid] = result['chain_unpID']
    
    # except Exception as e:
    #     print(f'{pdbid}: {e}')

# update multi_domain table
for no_unp in list_pdbChainID_no_unp:
    multi_domain = multi_domain[(multi_domain['FA-PDBID']!=no_unp[0])&(multi_domain['FA-CHAINID']!=no_unp[1])]
# save
save_tab(multi_domain, path_tab_scop_FA_continuous_uni_multi_parse_domian)

# save dict_pdb_chain_unp
dump_dict2json(dict_pdb_chain_unp, path_pdb_chain_unp)