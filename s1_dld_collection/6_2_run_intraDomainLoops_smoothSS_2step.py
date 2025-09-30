import os
import pandas as pd
import numpy as np
import json

from loop.dssp import save_chains, generate_pdp_chain_pairs
from utils.common import load_tab, dump_dict2json, int_domainREG, read_json2dict
from loop.extract_loop import check_terminus, check_overlap_2domain
from domain_linker.extract_linkers import get_2domain_start_end
from domain.scop2 import get_domain_start_end, get_domain_dataset

from params import *

'''
Extact Intra-domain Loops.
1. check if the loop overlapping with Inter-domain loops
2. check if the loop is N- ot C- terminus

3. save domain data: {'FA-DOMID': {}, ...}
    exclude terminus & exclude dependent/independent-domain-linkers
    start, end, ...
    based on 
        a. 1step-smoothing loops which are overlapped with two domains (no linker between 2-domains after 2step-smoothing.) - 0 loop like this.
        b. inter-domain-loops (2step-smoothing)
    to calculate/update
    
    Note: currently, only keep domains 
        a. in a multi-domain chain 
        b. with linkers (dependent- or independent-domain-linkers)
        c. with a 1step-smoothing loop which overlap with 2 domains (there is no linker between them). - 0 loop like this currently.
'''

# All loops that were extracted based on 2-SS smoothing
dict_all_loops_smoothSS_first = read_json2dict(path_loops_all_smoothSS_first)
# inter-domain loops
dict_domain_loops_smoothSS_second = read_json2dict(path_loops_domain_smoothSS_second)

# get Multi-domain table
multi_domain = load_tab(path_tab_scop_FA_continuous_uni_multi_parse_domian)
# get all domains' starts and ends for each dhain
# {(pdbid, chainid): [[domain1_start, domain1_end, domian1_ID], [domain2_start, domain2_end, domain2_ID], ...], }
dict_chain_domainID = get_domain_start_end(multi_domain)

# save Domain dataset: {'FA-DOMID': {'chain':, 'start':, 'end':, 'start_unp':, 'end_unp':, 'unp_acc':}}
dict_domain_dataset = get_domain_dataset(multi_domain)

num_intra_loops = 0

count_remove_linker = 0
count_remove_terminu_N = 0
count_remove_terminu_C = 0
count_remove_overlapOneDomain = 0
count_remove_overlapZeroDomain = 0
# loop through all chains
for pdb_chainid_unp, inter_loops in dict_domain_loops_smoothSS_second.items():
    # list of loops are going to be removed
    remove_loops = []
    # loop through all inter_domain_loops in a chain
    for interL in inter_loops:
        # loop through all intra_domain_loops
        for loop in dict_all_loops_smoothSS_first[pdb_chainid_unp]:
            # check if the current loop is inside an inter_domain loop
            # if Yes, delete the loop from dict_all_loops_smoothSS_2
            # if No, continue
            if (loop['start']>=interL['start']) & (loop['end']<=interL['end']):
                # save remove loops
                remove_loops.append(loop)
        '''
        preparing Domain Dataset:
        
        FA-DOMID: from scop
        {'FA-DOMID': {'pdb_chainid_unp':, 'chain':, 'start':, 'end':, 'start_unp':, 'end_unp':, 'unp_acc':}}
        '''
        # return {'d1_start', 'd1_end', 'd2_start', 'd2_end', 'd1_start_unp', 'd1_end_unp', 'd2_start_unp', 'd2_end_unp', 'd1_unp_acc', 'd2_unp_acc'}
        dict_2domain = get_2domain_start_end(interL, multi_domain)
        
        # Domain1
        if dict_domain_dataset[interL['domain1']]['end'] > dict_2domain['d1_end']:
            # update end & end_unp
            dict_domain_dataset[interL['domain1']]['end'] = dict_2domain['d1_end']
            dict_domain_dataset[interL['domain1']]['end_unp'] = dict_2domain['d1_end_unp']
            
        # Domain2
        if dict_domain_dataset[interL['domain2']]['start'] < dict_2domain['d2_start']:
            # update end & end_unp
            dict_domain_dataset[interL['domain2']]['start'] = dict_2domain['d2_start']
            dict_domain_dataset[interL['domain2']]['start_unp'] = dict_2domain['d2_start_unp']
           
    # remove loops
    count_remove_linker = count_remove_linker + len(remove_loops)
    for loop in remove_loops:
        dict_all_loops_smoothSS_first[pdb_chainid_unp].remove(loop)

# save N-/C-terminus
# format: {pdb_chainid_unp: [loops]}
dict_pdbChain_terminus = {}

# check if a loop is N- or C-terminus
for pdb_chainid_unp, loops in dict_all_loops_smoothSS_first.items():
    pdbid = pdb_chainid_unp.split('_')[0].upper()
    chainid = pdb_chainid_unp.split('_')[1].upper()
    unp = pdb_chainid_unp.split('_')[2]
    
    # save N-/C-terminus
    list_terminu = []
    
    # list of loops are going to be removed
    remove_loops = []
    if not ((pdbid, chainid) in dict_chain_domainID.keys()):
        print(f'chain does not exist: {(pdbid, chainid)}')
        continue
    list_domain = dict_chain_domainID[(pdbid, chainid)]
    for loop_idx in range(len(loops)):
        loop = loops[loop_idx]
        
        # check if this 'Intra-domain-loop' is overlapping with 2 domains.
        # {'overlap': , 'domains':[]}
        overlap_domain = check_overlap_2domain(loop, list_domain)
        
        terminu_domain = check_terminus(loop, list_domain)
        # 1. check temimus
        # 2. update domain dataset starts & ends
        
        if terminu_domain['terminu']:
            # 1. add to remove_list
            remove_loops.append(loop)
            # 2. update domain dataset
            # print(terminu_domain['domain'], terminu_domain['is-terminu'], terminu_domain['terminu'])
            domainData = dict_domain_dataset[terminu_domain['domain']]
            
            # N-terminu
            if terminu_domain['tag']=='N':
                domainData['start'] = loop['end'] + 1
                domainData['start_unp'] = loop['end_unp'] + 1
                count_remove_terminu_N = count_remove_terminu_N + 1
                # save terminus
                loop['terminu'] = 'N'
                list_terminu.append(loop)
            # C-terminu
            else:
                domainData['end'] = loop['start'] - 1
                domainData['end_unp'] = loop['start_unp'] - 1
                count_remove_terminu_C = count_remove_terminu_C + 1
                
                # save terminus
                loop['terminu'] = 'C'
                list_terminu.append(loop)
                
            dict_domain_dataset[terminu_domain['domain']] = domainData
        # intra-domain loops & loops overlap with 0 or 1 domain
        else:
            if terminu_domain['tag']=='intra':
                # update loop with domain info
                loop['domain'] = terminu_domain['domain']
                dict_all_loops_smoothSS_first[pdb_chainid_unp][loop_idx] = loop
            elif terminu_domain['tag']=='overlap1_start':
                domainData = dict_domain_dataset[terminu_domain['domain']]
                domainData['start'] = loop['end'] + 1
                domainData['start_unp'] = loop['end_unp'] + 1
                
                dict_domain_dataset[terminu_domain['domain']] = domainData
                # remove
                remove_loops.append(loop)
                count_remove_overlapOneDomain = count_remove_overlapOneDomain + 1
            elif terminu_domain['tag']=='overlap1_end':
                domainData = dict_domain_dataset[terminu_domain['domain']]
                domainData['end'] = loop['start'] - 1
                domainData['end_unp'] = loop['start_unp'] - 1
                
                dict_domain_dataset[terminu_domain['domain']] = domainData
                # remove
                remove_loops.append(loop)
                count_remove_overlapOneDomain = count_remove_overlapOneDomain + 1
            else: 
                remove_loops.append(loop)
                count_remove_overlapZeroDomain = count_remove_overlapZeroDomain + 1
    # save terminus
    if len(list_terminu)>0:
        dict_pdbChain_terminus[pdb_chainid_unp] = list_terminu
    
    # remove loops            
    for loop in remove_loops:
        dict_all_loops_smoothSS_first[pdb_chainid_unp].remove(loop)

print(f'count_remove_linker: {count_remove_linker}')
print(f'count_remove_terminu_N: {count_remove_terminu_N}')
print(f'count_remove_terminu_C: {count_remove_terminu_C}')
print(f'count_remove_overlapOneDomain: {count_remove_overlapOneDomain}')
print(f'count_remove_overlapZeroDomain: {count_remove_overlapZeroDomain}')

num_intra_loops = sum([len(loops) for pdb_chainid_unp, loops in dict_all_loops_smoothSS_first.items()])
print(f'The number of intra-domain loops: {num_intra_loops}')
print('Saving intra-domain loops ...')
dump_dict2json(dict_all_loops_smoothSS_first, path_loops_intra_domain_smoothSS_first_2step)

print('Saving N-/C-terminus ...')
dump_dict2json(dict_pdbChain_terminus, path_terminus)

print('Save domain dataset ...')
dump_dict2json(dict_domain_dataset, path_domain_dataset_pre)

print('Done!')