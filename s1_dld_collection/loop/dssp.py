import os
import pandas as pd
import numpy as np

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBList

from utils.protein import protein_letters_3to1
from utils.common import dump_dict2json
from loop.pdb_files import get_dssp_key

###
# DSSP: aim to find loops in multi-domain proteins
# tool: Biopython.DSSP
# date: 11 Mar, 2023
# Functions related to extract loops by using DSSP
# !!!Note: even if mmCif is the current achive PDB file, but in this context, .pdb files are used. This is beacuse that DSSP deal with .pdb file.
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

def map_8ssTo3ss(SS8_seq: list) -> list:
    '''
    Convert DSSP's 8-state assignments into 3-state [C - coil, E - extended (beta-strand), H - helix
        The DSSP codes for 8-state secondary structure used here are:
        # =====     ====
        # Code      Structure
        # =====     ====
        # H         Alpha helix (4-12)
        # B         Isolated beta-bridge residue
        # E         Strand
        # G         3-10 helix
        # I         Pi helix
        # T         Turn
        # S         Bend
        # -         None
        # =====     ====
        
    params:
        SS8_seq - str, 8-state secondary structure
        
    return:
        SS8_seq - str, 3-state secondary structure
    '''
    
    '''
    SS-Scheme 1: H,G,I->H ; E,B->E ; T,S->C    The one we use here, in paper5_00
    SS-Scheme 2: H,G->H ; E,B->E ; I,T,S->C    I think this is most common
    
    SS-Scheme 3: H,G->H ; E->E ; I,B,T,S->C
    
    SS-Scheme 4: H->H ; E,B->E ; G,I,T,S->C
    SS-Scheme 5: H->H ; E->E ; G,I,B,T,S->C
    '''
    
    '''
    SS-Scheme we use: 
        H,G,I->H ; E->E; B->B ; _,T,S->C
    '''
    
    dict_map = {'-': 'C', 'T': 'C', 'S': 'C', 
                'G': 'H', 'H': 'H', 'I': 'H', 
               'E': 'E', 'B': 'B'}
    # SS3_seq = ''.join([dict_map[x] for x in SS8_seq])
    SS3_seq = [dict_map[x] for x in SS8_seq]
    return SS3_seq
    
def save_SS(pdbIDs: list, dir_pdb: str, path_dssp_pdbid: str):
    '''
    Using DSSP, get secondary structure from pdb files. Save information to individual json file (per PDBID)
    https://biopython.org/docs/1.76/api/Bio.PDB.DSSP.html
    
    params:
        pdbIDs - list, a list of pdbids.
        dir_pdb - directory to PDB files.
        path_dssp_pdbid - directory to save pdb_secondary_structure info.
        
    return:
        error_pids - list of pdbids that cannot be parsed by DSSP.
    Note:
        ss - {pdbid: '', 'sequence':, 'chain': , 'pdb_seqidx':, 
                                        'sec_structure_8': , 'sec_structure_3':}
    '''
    
    '''
    dssp file format:
    
    Tuple Index - Value
    0 - DSSP index
    1 - Amino acid
    2 - Secondary structure
    3 - Relative ASA
    4 - Phi
    5 - Psi
    6 - NH–>O_1_relidx
    7 - NH–>O_1_energy
    8 - O–>NH_1_relidx
    9 - O–>NH_1_energy
    10 - NH–>O_2_relidx
    11 - NH–>O_2_energy
    12 - O–>NH_2_relidx
    13 - O–>NH_2_energy
    '''
    
    # pdbids that cannot be parsed by DSSP
    error_pids = []
    
    # parse structure
    cifP = MMCIFParser()
    
    # loop through all pdbids, get all secondary structures
    for pid in pdbIDs:
        pid = pid.lower()
        
        try: 
            structure = cifP.get_structure(pid, dir_pdb+'/%s.cif' % pid)

            # use only the first model
            model = structure[0]

            # calculate DSSP
            # file_type also could be 'PDB'
        
            dssp = DSSP(model, dir_pdb+'/%s.cif' % pid, file_type='mmCif')

            # extract sequence and secondary structure from the DSSP tuple
            sequence = []
            sec_structure_8 = []
            sec_structure_3 = []
            chain = []
            pdb_seqidx = []

            # dssp index, for locating Hbond
            dssp_idx = []
            
            # ('D', (' ', 188, ' '))
            dssp_key = []
            
            # H-bond
            NH2O_1_relidx = []
            NH2O_1_energy = []
            O2NH_1_relidx = []
            O2NH_1_energy = []
            NH2O_2_relidx = []
            NH2O_2_energy = []
            O2NH_2_relidx = []
            O2NH_2_energy = []

            for i in range(len(dssp)):
                # https://bioinformatics.stackexchange.com/questions/11587/what-is-the-aim-of-insertion-codes-in-the-pdb-file-format
                # ignore the insertion residues
                if dssp.property_keys[i][1][2]!=' ':
                    continue
                a_key = list(dssp.keys())[i]

                sequence.append(dssp[a_key][1])
                sec_structure_8.append(dssp[a_key][2])

                # dssp idx
                dssp_idx.append(dssp[a_key][0])
                
                # ('D', (' ', 188, ' '))
                dssp_key.append(a_key)
                
                # H-bond
                NH2O_1_relidx.append(dssp[a_key][6])
                NH2O_1_energy.append(dssp[a_key][7])
                O2NH_1_relidx.append(dssp[a_key][8])
                O2NH_1_energy.append(dssp[a_key][9])
                NH2O_2_relidx.append(dssp[a_key][10])
                NH2O_2_energy.append(dssp[a_key][11])
                O2NH_2_relidx.append(dssp[a_key][12])
                O2NH_2_energy.append(dssp[a_key][13])

                chain.append(dssp.property_keys[i][0])
                pdb_seqidx.append(str(dssp.property_keys[i][1][1]))

            # convert DSSP's 8-state assignments into 3-state [C - coil, E - extended (beta-strand), H - helix]
            sec_structure_3 = map_8ssTo3ss(sec_structure_8)

            # save values to dict
            ss = {'dssp_key': dssp_key, 'dssp_idx': dssp_idx, 'pdbid': pid, 'sequence': sequence, 'chain': chain, 
                  'auth_seqidx': pdb_seqidx, 'sec_structure_8': sec_structure_8, 'sec_structure_3': sec_structure_3, 
                       'NH2O_1_relidx': NH2O_1_relidx, 'NH2O_1_energy': NH2O_1_energy,
                      'O2NH_1_relidx': O2NH_1_relidx, 'O2NH_1_energy': O2NH_1_energy,
                      'NH2O_2_relidx': NH2O_2_relidx, 'NH2O_2_energy': NH2O_2_energy,
                      'O2NH_2_relidx': O2NH_2_relidx, 'O2NH_2_energy': O2NH_2_energy,}
            
            # save dict to [pdbid].json
            file_name = pid+'.json'
            path_json = os.path.join(path_dssp_pdbid, file_name)
            dump_dict2json(ss, path_json)
            
        except Exception as e:
            print(f'{pid}: {e}')
            error_pids.append([pid, e])
    
    return error_pids

def save_chains_full(pdbid: str, df_ss: pd.DataFrame, dir_chain: str, dict_pdb_chain: dict, path_pdbe_uniprot_pdbid:str):
    '''
    Group chains of a PDB entry. save each chain into an JSON file.
    
    Params:
        pdbid - e.g. '1cdl'
        df_ss - information of an PDB entry in a table
        dir_chain - folder to save JSON files, file format: [pdbid]_[chainid].json
        dict_pdb_chain - {pdbid: [chainids]}
        path_pdbe_uniprot_pdbid - df to json file path, pdb entry information from pdbe, including Uniprot information.
        
    return:
        result - {'list_pdbChainID_no_unp': [], 'pdbid': {chainid: [unpID]}}
            pdb_chain_unp - {chainid: [unpID]}
            list_pdbChainID_no_unp - list of [[pdbid, chainid], ....] cannot map to Uniprot.
    '''
    # print(pdbid)
    pdbid = pdbid.lower()
    df_pdbe_uniprot = pd.read_json(os.path.join(path_pdbe_uniprot_pdbid, pdbid+'.json'))
    df_ss_group = df_ss.groupby('chain')
    # group pdbe data based on chainID
    df_pdbe_uniprot_group = df_pdbe_uniprot.groupby('chain')
    # uniprot_columns = ['asym_id', 'seq_id', 'unp_res', 'unp_num', 'unp_acc', 'unp_segment_id', 'unp_instance_id']
    # save chains cannot map to Uniprot
    list_pdbChainID_no_unp = []
    chain_unpID = {}
    
    for chainid in dict_pdb_chain[pdbid.upper()]:
        try:
            if chainid in df_ss_group.groups.keys():
                df_chain = df_ss_group.get_group(chainid)
                df_chain_uniprot = df_pdbe_uniprot_group.get_group(chainid)
            else:
                df_chain = df_ss_group.get_group(int(chainid))
                df_chain_uniprot = df_pdbe_uniprot_group.get_group(int(chainid))
        except:
            list_pdbChainID_no_unp.append([pdbid, chainid])
            continue
            
        ''' 
        if df_chain.shape[0]!=df_chain_uniprot.shape[0]:
            print(f'size is not equal: {pdbid, chainid}')
        '''
        # merge
        df_concat = df_chain_uniprot.merge(df_chain, on=['dssp_key_str', 'chain', 'auth_seqidx'], how='left', validate='one_to_one')
        
        # fill na, some residues are not in DSSP
        df_concat['dssp_idx'] = df_concat['dssp_idx'].fillna(-1)
        df_concat['pdbid'] = df_concat['pdbid'].fillna(method='ffill')
        
        df_concat['sequence'] = df_concat['sequence'].fillna(df_concat['unp_res'])
        df_concat['sec_structure_8'] = df_concat['sec_structure_8'].fillna('-')
        df_concat['sec_structure_3'] = df_concat['sec_structure_3'].fillna('C')
        
        df_concat['missing'] = df_concat['missing'].fillna(0)
        
        df_concat = df_concat.fillna(0)
        
        df_concat_group = df_concat.groupby('unp_acc')
        
        for unpID in df_concat_group.groups.keys():
            
            df_unp = df_concat_group.get_group(unpID)
            # file name format: [pdbid]_[chainid].json
            chain_path = os.path.join(dir_chain, f'{pdbid}_{chainid.lower()}_{unpID}.json')
            df_unp.to_json(chain_path)
        
        # save UNIPROT ID
        chain_unpID[chainid] = list(df_concat_group.groups.keys())
        
        # # one chain map to more than one UniprotID
        # if len(df_concat['unp_acc'].unique())>1:
        #     list_pdbChainID_no_unp.append([pdbid, chainid])
    return {'list_pdbChainID_no_unp': list_pdbChainID_no_unp, 'chain_unpID': chain_unpID}

def save_chains(pdbid: str, df_ss: pd.DataFrame, dir_chain: str, dict_pdb_chain: dict):
    '''
    Group chains of a PDB entry. save each chain into an JSON file.
    
    Params:
        pdbid - e.g. '1cdl'
        df_ss - information of an PDB entry in a table
        dir_chain - folder to save JSON files, file format: [pdbid]_[chainid].json
        dict_pdb_chain - {pdbid: [chainids]}
    '''
    pdbid = pdbid.lower()
    
    df_ss_group = df_ss.groupby('chain')
    
    for chainid in dict_pdb_chain[pdbid.upper()]:
        try:
            df_chain = df_ss.groupby('chain').get_group(chainid)
        except:
            df_chain = df_ss.groupby('chain').get_group(int(chainid))
        
        # file name format: [pdbid]_[chainid].json
        chain_path = os.path.join(dir_chain, pdbid+'_'+chainid.lower()+'.json')
        df_chain.to_json(chain_path)
        
def generate_pdp_chain_pairs(df_domain: pd.DataFrame):
    '''
    get {pdbid: [chains]} pair.
    
    params:
        df_domain - scop domain table, generate from SCOP-related code
        
    return: 
        dict_pdb_chain - {pdbid: [chainids]}
    '''
    
    list_pdbid = df_domain['FA-PDBID'].unique()
    
    group_domain = df_domain.groupby(['FA-PDBID'])
    dict_pdb_chain = {}

    for pdbid in list_pdbid:
        chains = list(group_domain.get_group(pdbid)['FA-CHAINID'].unique())
        dict_pdb_chain[pdbid] = chains
    
    return dict_pdb_chain

def add_missing_residues(pdbid: str, df_ss: pd.DataFrame, df_missing: pd.DataFrame) -> pd.DataFrame:
    '''
    comlete the dssp/secondary structure with adding missing residues.
    Add dssp_key_str to pdb entries.
    
    params:
        df_ss - table without missing residues
        df_missing - pd.DataFrame, missing residues get from mmcif file
                     if it is empty, no missing residue
    return:
        df_all - table with missing residues.
    '''
    df_ss['missing'] = 0
    
    # no missing rresidue
    if df_missing.empty:
        print(f'{pdbid}: No missing residues')
        df_ss['dssp_key_str'] = [str(key[0])+'_'+str(key[1][1])+'_'+str(key[1][2]) for key in df_ss['dssp_key']]
        return df_ss
    
    df_missing['dssp_idx'] = -1
    df_missing['pdbid'] = pdbid
    df_missing['auth_comp_id'] = protein_letters_3to1(df_missing['auth_comp_id'])

    df_missing['sec_structure_8'] = '-'
    df_missing['sec_structure_3'] = 'C'
    df_missing['missing'] = 1

    # no H-bond information
    # [-1, -40]
    df_missing['NH2O_1_relidx'] = 0
    df_missing['NH2O_1_energy'] = 0
    df_missing['O2NH_1_relidx'] = 0
    df_missing['O2NH_1_energy'] = 0
    df_missing['NH2O_2_relidx'] = 0
    df_missing['NH2O_2_energy'] = 0
    df_missing['O2NH_2_relidx'] = 0
    df_missing['O2NH_2_energy'] = 0
    
    df_missing = df_missing.rename(columns={'auth_asym_id': 'chain', 'auth_comp_id': 'sequence', 'auth_seq_id': 'auth_seqidx'})
    
    # remove those residues which are not missing residue from df_missing
    df_missing["auth_seqidx"] = pd.to_numeric(df_missing["auth_seqidx"])

    # get existing chain,index pair
    df_missing['dssp_key_str'] = [str(key[0])+'_'+str(key[1][1])+'_'+str(key[1][2]) for key in df_missing['dssp_key']]
    df_ss['dssp_key_str'] = [str(key[0])+'_'+str(key[1][1])+'_'+str(key[1][2]) for key in df_ss['dssp_key']]
    
    df_missing = df_missing[~df_missing['dssp_key_str'].isin(df_ss['dssp_key_str'])]
    print(pdbid, df_missing.shape[0])
    # convert column "a" of a DataFrame
    df_all = pd.concat([df_ss, df_missing])
    df_all["auth_seqidx"] = pd.to_numeric(df_all["auth_seqidx"])

    df_all = df_all.sort_values(['chain', 'auth_seqidx']).reset_index(drop=True)
    
    return df_all