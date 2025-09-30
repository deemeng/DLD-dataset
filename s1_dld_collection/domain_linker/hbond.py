import os
import pandas as pd
import numpy as np

from params import *

'''
Functions to find and record H-bonds.
# !!!!bond between a-b and b-a is the same!!!!

Di
20, Mar, 2023
'''

def check_2Hbond(df_chain: pd.DataFrame, loop: dict, dict_2domains:dict):
    '''
    Check if there are 2 or more H-bonds between 2 domains
    consider a stable connection between 2 domains if there are 2 or more H-honds.
    
    params:
        df_chain - chain infomation generated from DSSP file
        loop - domain-loop information
        dict_2domains - {d1_start:, d1_end, d2_start, d2_end}, domain1 and domain2 start&end information
    return:
        hbond_idx - list of hbond indices
    '''
    
    list_energy = ['NH2O_1_energy', 'O2NH_1_energy', 'NH2O_2_energy', 'O2NH_2_energy']
    list_relidx = ['NH2O_1_relidx', 'O2NH_1_relidx', 'NH2O_2_relidx', 'O2NH_2_relidx']
    
    # consider a stable connection between 2 domains if there are 2 or more H-honds    
    # save H-bonds info
    hbond_idx = set()
    for relidx, energy in zip(list_relidx, list_energy):
        # df_hbond = df_chain[df_chain[energy].between(-40, -1, inclusive='both')]
        df_hbond = df_chain[df_chain[energy]<-1]
        
        for _, row in df_hbond.iterrows():
            authIdx = int(row['auth_seqidx'])
            # real
            rel_index = int(row[relidx]) + int(row['dssp_idx'])
            
            try:
                rel_authIdx = int(df_chain.loc[df_chain['dssp_idx']==rel_index, 'auth_seqidx'].item())
            except:
                # connection between the current chain and another chain
                continue

            # check the bond location between two domains or not
            if (min(authIdx, rel_authIdx)<loop['start']) & (min(authIdx, rel_authIdx)>=dict_2domains['d1_start']) \
                & (max(authIdx, rel_authIdx)>loop['end']) & (max(authIdx, rel_authIdx)<=dict_2domains['d2_end']):
                # !!!!bond between a-b and b-a is the same!!!!
                print('relative idx')
                print((relidx, rel_index))
                print((authIdx, rel_authIdx))
                hbond_idx.add((authIdx, rel_authIdx))
                hbond_idx.add((rel_authIdx, authIdx))
                
    return list(hbond_idx)