import os
import numpy as np
import pandas as pd

import json

from params import *
from utils.common import load_tab, df2json, save_tab, read_json2dict, dump_dict2json

'''
complete information (columns) for domain data.

{
domain1: 
    {
    pdb_chainid, unp_acc,
    start, end, start_unp, end_unp, length,
    miss_length, miss_percentage, missing,
    seq, seq_unp, seq_id, seq_id_unp, dssp_key_str
    },
domain2: {}, ..
}

Di
7 May 2023
'''
# save keys cannot map to Uniprot.
remove_list = []

multi_domain = load_tab(path_tab_scop_FA_continuous_uni_multi_parse_domian)
dict_domain_dataset = read_json2dict(path_domain_dataset_pre)

for key, domain in dict_domain_dataset.items():
    pid_chain = domain['pdb_chainid']
    unp_acc = domain['unp_acc']
    try:
        df_chain = pd.read_json(os.path.join(path_dssp_full_uniprot_chainid, f'{pid_chain}_{unp_acc}.json'))
    except:
        remove_list.append(key)
        print(f'{key}: cannot map to Uniprot')
        continue
    df_domain = df_chain[(df_chain['unp_acc']==domain['unp_acc']) & (df_chain['unp_num']>=domain['start_unp']) & (df_chain['unp_num']<=domain['end_unp'])]
    if df_domain.shape[0]==0:
        remove_list.append(key)
        print(f'{key}: cannot map to Uniprot')
        continue
    
    domain['dssp_key_str'] = list(df_domain['dssp_key_str'])
    domain['length'] = df_domain.shape[0]
    domain['miss_length'] = sum(df_domain['missing'])
    domain['miss_percentage'] = sum(df_domain['missing'])/df_domain.shape[0]
    domain['missing'] = list(df_domain['missing'])
    domain['seq'] = list(df_domain['sequence'])
    domain['seq_unp'] = list(df_domain['unp_res'])
    domain['seq_id'] = list(df_domain['auth_seqidx'])
    domain['seq_id_unp'] = list(df_domain['unp_num'])
    
    domain['start'] = list(df_domain['auth_seqidx'])[0]
    domain['end'] = list(df_domain['auth_seqidx'])[-1]
    domain['start_unp'] = list(df_domain['unp_num'])[0]
    domain['end_unp'] = list(df_domain['unp_num'])[-1]
    
    dict_domain_dataset[key] = domain
    
# remove items/domains cannot map to Uniprot.
for key in remove_list:
    dict_domain_dataset.pop(key)
    
print('save domain dataset ...')
dump_dict2json(dict_domain_dataset, path_domain_dataset)