import os
import numpy as np
import pandas as pd

import json
from pdbecif.mmcif_tools import MMCIF2Dict as MMCIF2Dict_pdbe

from params import *

from utils.common import load_tab, df2json, save_tab
from loop.dssp import generate_pdp_chain_pairs
from loop.pdb_files import get_uniprot_info
'''
25 Apr, 2023
Di

Aim to add Uniprot infomation (from pdbe sifts_columns) to path_dssp_full_pdbid
and save the results to folder path_dssp_full_uniprot_pdbid
'''
multi_domain = load_tab(path_tab_scop_FA_continuous_uni_multi_parse_domian)
dict_pdb_chain = generate_pdp_chain_pairs(multi_domain)

no_uniprot_mapping = []

for pdbid in list(dict_pdb_chain.keys()):
    # print(pdbid)
    dict_mmcif_pdbe = MMCIF2Dict_pdbe()
    dict_cif_pdbe = dict_mmcif_pdbe.parse(os.path.join(path_pdbe_cif, f'{pdbid.lower()}_updated.cif'))
    
    # get uniprot info
    df_uniprot = get_uniprot_info(pdbid, dict_cif_pdbe)
    
    # check if 
    if df_uniprot.shape[0]==0:
        no_uniprot_mapping.append(pdbid)
        continue
    
    # merge
    try:
        # for numeric chainID 1, 2, ... instead of letters
        df_uniprot = df_uniprot.astype({'chain': int})
        df_uniprot = df_uniprot.sort_values(['chain', 'auth_seqidx']).reset_index(drop=True)
    except:
        df_uniprot = df_uniprot.sort_values(['chain', 'auth_seqidx']).reset_index(drop=True)
        
    # save df_all to file
    df2json(df_uniprot, os.path.join(path_pdbe_uniprot_pdbid, pdbid.lower()+'.json'))

print(f'PDB entries cannot map to Uniprot: {no_uniprot_mapping}')
print('Update path_tab_scop_FA_continuous_uni_multi_parse_domian')
# ['6S9U']
multi_domain = multi_domain[~multi_domain['FA-PDBID'].isin(no_uniprot_mapping)]
# save
save_tab(multi_domain, path_tab_scop_FA_continuous_uni_multi_parse_domian)
print('Done!')