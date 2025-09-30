import os
import pandas as pd
import numpy as np

from loop.pdb_files import get_PDBfiles, parse_mmcif
from loop.dssp import add_missing_residues
from utils.common import load_tab, save_tab, dump_dicts2jsons, lower_pdbid, df2json
from utils.protein import protein_letters_3to1
from loop.pdb_files import get_dssp_key
from params import *

'''
17 Mar, 2023
Di

Aim to complete dssp/pdbid file with missing residues

Notice:
 1. in pdb file '_pdbx_unobs_or_zero_occ_residues' includes not onlyt the missing residues but also some others. 
 2. Therefore, check if the residue is already in the chains before adding missing residues.
'''
# get PDBIDs
multi_domain = load_tab(path_tab_scop_FA_continuous_uni_multi_parse_domian)
list_PDBID = multi_domain['FA-PDBID'].unique()
list_pdbid = lower_pdbid(list_PDBID)

# protein entries without missing residues
list_nomissing = []

print(f'adding missing resiues, and saving to {path_dssp_full_pdbid}')
# loop through all .json file and add missing residues
# save them to dssp/dull_pdbid
for pdbid in list_pdbid:
    ss_path = os.path.join(path_dssp_pdbid, pdbid + '.json')
    df_ss = pd.read_json(ss_path)
    # get missing residues
    # information needed from mmcif file
    path_cif = os.path.join(path_pdb_cif, pdbid + '.cif')
    pre_key = '_pdbx_unobs_or_zero_occ_residues.'
    tail_keys = ['PDB_model_num', 'auth_asym_id', 'auth_comp_id', 'auth_seq_id', 'PDB_ins_code']
    
    # with or without missing residues
    try:
        df_missing = parse_mmcif(path_cif, pre_key, tail_keys)
        # NMR experiments have more than one models
        # Only take the 1st one
        df_missing = df_missing[df_missing['PDB_model_num']=='1']
        list_ins_code = [' ' if ins_code=='?' else ins_code for ins_code in df_missing['PDB_ins_code']]
        dssp_key = get_dssp_key(df_missing['auth_asym_id'], df_missing['auth_seq_id'], list_ins_code)
        df_missing['dssp_key'] = dssp_key
        df_missing.drop(['PDB_model_num', 'PDB_ins_code'], inplace=True, axis=1)
    except:
        list_nomissing.append(pdbid)
        df_missing = pd.DataFrame()
        
    df_all = add_missing_residues(pdbid, df_ss, df_missing)
    
    # save df_all to file
    df2json(df_all, os.path.join(path_dssp_full_pdbid, pdbid+'.json'))

print('Save pdb entries with completing missing residues...')
# save pdbids without missing residue
with open(path_nomissing_residue, 'w') as f:
    for line in list_nomissing:
        f.write(f"{line}\n")
print('Done!')