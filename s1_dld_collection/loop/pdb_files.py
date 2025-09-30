import os
import pandas as pd
import numpy as np

import wget

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBList
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
###
# DSSP: aim to find loops in multi-domain proteins
# tool: Biopython.DSSP
# date: 11 Mar, 2023
# Functions related to extract loops by using DSSP
# !!!Note: even if mmCif is the current achive PDB file, but in this context, .pdb files are used. This is beacuse that DSSP deal with .pdb file.
###

###
# RCSB PDB
###
def get_PDBfiles(pdbIDs: list, out_dir: str):
    '''
    params:
        pdbIDs - list, a list of pdbids.
        out_dir - directory to save PDB files
        
    return:
        None
    '''
    
    # initialize PDB downloader
    pdb_dl = PDBList()
    
    # download PDB files, save them in current directory
    for pid in pdbIDs:
        # https://biopython.org/docs/1.76/api/Bio.PDB.PDBList.html
        # file_format='pdb'
        pdb_dl.retrieve_pdb_file(pid, pdir=out_dir, file_format='mmCif', overwrite=True)
        
def parse_mmcif(path_cif, pre_key: str, tail_keys: list) -> pd.DataFrame: 
    '''
    Get infomation from mmcif file according to the keys
    
    params:
        path_cif - path to .cif file
        pre_key - str, e.g. _pdbx_unobs_or_zero_occ_residues.
        tail_keys - list, e.g. ['auth_asym_id', 'auth_comp_id', 'auth_seq_id']. pre_key+tail_key is the keys/columns
    return:
        df_all - table with keys as columns.
    '''
    
    dict_cif  = MMCIF2Dict(path_cif)
    df_cif = pd.DataFrame({k: dict_cif.get(pre_key+k) for k in tail_keys})
    
    return df_cif

###
# PDBe
###
def download_pdbe_updatedCif(list_pdbid: list, path_pdbe_cif):
    '''
    Download updated mmcif files from PDBe, which includes SIFTS infomation (uniprot ID)
    
    params:
        list_pdbid - list of pdbIDs need to download. 
        path_pabe_cif - folder to save those mmcif files
    '''
    for pdbid in list_pdbid:
        URL = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdbid.lower()}_updated.cif"
        response = wget.download(URL, path_pdbe_cif)

def get_dssp_key(l_chain: list, l_auth_seqidx: list, l_ins_code: list):
    '''
    Given chainID, auth_seqidx, and PDB_ins_code, generate DSSP key list.
    DSSP key: e.g. ('B', (' ', 215, ' '))
    
    params:
        l_chain - list of chainIDs (author)
        l_auth_seqidx - list of auth_seqidx
        l_ins_code - list of PDB_ins_code
        
    return:
        dssp_key - [('B', (' ', 215, ' ')), ....]
    '''
    dssp_key = [(chainid, (' ', auth_seqidx, ins_code)) for (chainid, auth_seqidx, ins_code) in 
     zip(l_chain, l_auth_seqidx, l_ins_code)]
    
    return dssp_key
    
def get_uniprot_info(pdbid: str, dict_cif_pdbe: dict) -> pd.DataFrame:
    '''
    Given a dict from pdbecif.mmcif_tools.MMCIF2Dict object, get Uniprot information.
    
    params:
        pdbid - pdb entry ID.
        dict_cif_pdbe - dict, generate from MMCIF2Dict_pdbe.parse()
    return:
        df_uniprot - including columns: 'asym_id', 'seq_id', 'unp_res', 'unp_num', 'unp_acc', 'unp_segment_id',
       'unp_instance_id', 'auth_seqidx', 'chain'
    '''
    
    ###
    # 1. Uniprot information from SIFTS columns
    ###
    try:
        sifts_cols = ['asym_id', 'seq_id', 'unp_res', 'unp_num', 'unp_acc', 'unp_segment_id', 'unp_instance_id']
        # cif_dict
        df_pdbe_sift = pd.DataFrame(dict_cif_pdbe[pdbid]['_pdbx_sifts_xref_db'])[sifts_cols]
        df_pdbe_sift = df_pdbe_sift.drop_duplicates()
        
    except:
        return pd.DataFrame()
    
    ###
    # 2. get extra columns chain & auth_seqidx. For future using to connect with pdb_entry
    ###
    
    info_cols = ['asym_id', 'seq_id', 'pdb_seq_num', 'pdb_strand_id', 'pdb_ins_code']
    df_connect_info = pd.DataFrame(dict_cif_pdbe[pdbid.upper()]['_pdbx_poly_seq_scheme'])[info_cols]
    list_pdb_ins_code = [' ' if ins_code=='.' else ins_code for ins_code in df_connect_info['pdb_ins_code']]
    # DSSP key: ('B', (' ', 215, ' '))
    dssp_key = get_dssp_key(df_connect_info['pdb_strand_id'], df_connect_info['pdb_seq_num'], list_pdb_ins_code)
    
    df_connect_info['dssp_key'] = dssp_key
    df_connect_info = df_connect_info.drop('pdb_ins_code', axis=1)
    ###
    # 3. merge UniprotInfo & connectInfo
    ###
    df_uniprot = df_connect_info.merge(df_pdbe_sift, on=['asym_id', 'seq_id'], how='left')
    df_uniprot = df_uniprot.rename(columns={'pdb_strand_id': 'chain', 'pdb_seq_num': 'auth_seqidx'})
    # change datatype
    df_uniprot = df_uniprot.replace('?', np.nan)
    df_uniprot = df_uniprot.astype({'seq_id': int, 'unp_num': float, 'auth_seqidx': int})
    df_uniprot = df_uniprot.dropna()
    df_uniprot = df_uniprot.drop_duplicates()
    
    df_uniprot['dssp_key_str'] = [str(key[0])+'_'+str(key[1][1])+'_'+str(key[1][2]) for key in df_uniprot['dssp_key']]
    
    return df_uniprot