import os
import pandas as pd
import numpy as np

from utils.common import int_domainREG, dump_dict2json



def get_domain_loops(multi_domain: pd.DataFrame, loops: dict):
    '''
    Finds all loops that overlap with domains.
    
    params:
        multi_domain - table of all domain info, from SCOP
        loops - dict, all loops info, {'pid_chainid':{'start', 'end', 'seq'}}
        
    return:
        domain_loops - dict, all domain loops. {'pid_chainid':{'start', 'end', 'seq', domain1, domain2}}
        adjacent_domains - list, all domains collected by domain-loops
        num_loops - total number of loops.
    '''
    
    # col_names = ['FA-PDBREG-START', 'FA-PDBREG-END', 'FA-UNIREG-START', 'FA-UNIREG-END']
    # col_names = ['FA-PDBREG-START', 'FA-PDBREG-END']
    # multi_domain = int_domainREG(multi_domain, col_names)
    
    multi_domain = multi_domain.sort_values(['FA-PDBID', 'FA-CHAINID', 'FA-PDBREG-START'])
    
    group_chain = multi_domain.groupby(['FA-PDBID', 'FA-CHAINID'])
    dict_groups = group_chain.groups
    
    # all domain loops
    domain_loops = {}

    adjacent_domains = []

    # total number of loops
    num_loops = 0
    
    # UNIPROT does not match
    remove_pdbchain = []
    for key, idx in dict_groups.items():

        # domain_loop for a chain
        chain_domain_loops = []

        for i in range(1, len(idx)):
            # Linker connect 2 domains
            # the first domain
            domain1_start = int(multi_domain.loc[idx[i-1], 'FA-PDBREG-START'])
            domain1_end = int(multi_domain.loc[idx[i-1], 'FA-PDBREG-END'])
            # the second domain
            domain2_start = int(multi_domain.loc[idx[i], 'FA-PDBREG-START'])
            domain2_end = int(multi_domain.loc[idx[i], 'FA-PDBREG-END'])
            # uniprot_acc
            domain1_unp = str(multi_domain.loc[idx[i-1], 'FA-UNIID'])
            domain2_unp = str(multi_domain.loc[idx[i], 'FA-UNIID'])
            
            if domain1_unp!=domain2_unp:
                continue
                
            pid_chain = key[0].lower()+'_'+key[1].lower()+'_'+domain1_unp
            if not (pid_chain in loops.keys()):
                print(f'key - {pid_chain}: does no match')
                
                remove_pdbchain.append([key[0].upper(), key[1].upper()])
                continue
                
            loops4chain = loops[pid_chain]
            for loop in loops4chain:
                # overlap with 2 adjacent domains
                if (loop['start']<=domain1_end) & (loop['start']>domain1_start) & (loop['end']>=domain2_start) & (loop['end']<domain2_end):
                    # one domain loop
                    d_loop = loop
                    d_loop['domain1'] = int(multi_domain.loc[idx[i-1], 'FA-DOMID'])
                    d_loop['domain2'] = int(multi_domain.loc[idx[i], 'FA-DOMID'])

                    chain_domain_loops.append(d_loop)

                    adjacent_domains.append(int(multi_domain.loc[idx[i-1], 'FA-DOMID']))
                    adjacent_domains.append(int(multi_domain.loc[idx[i], 'FA-DOMID']))

                    num_loops = num_loops + 1
        if len(chain_domain_loops) > 0:            
            domain_loops[pid_chain] = chain_domain_loops
            
    return domain_loops, list(set(adjacent_domains)), remove_pdbchain, num_loops


def get_2domain_start_end(loop: dict, multi_domain: pd.DataFrame) -> dict:
    '''
    Get 2 adjacent domains' starts and ends, according the loop info.
    Domain1_start & Domain2_end from scoop.
    Domain1_end = loop_start - 1
    Domain2_start = loop_end + 1.
    
    params:
        loop - dict, inter-domain Loops, {'start', 'end', 'seq', 'domain1_domainID', 'domain2_domainID'}
        multi_domain - domain information from scoop
    return:
        dict_2domains - {'d1_start', 'd1_end', 'd2_start', 'd2_end', 'd1_start_unp', 'd1_end_unp', 'd2_start_unp', 'd2_end_unp', 
                        'd1_unp_acc', 'd2_unp_acc'}
    '''
    
    dict_2domains = {}
    
    domain1 = loop['domain1']
    domain2 = loop['domain2']
    
    dict_2domains['d1_unp_acc'] = multi_domain.loc[multi_domain['FA-DOMID']==domain1, 'FA-UNIID'].item()
    dict_2domains['d2_unp_acc'] = multi_domain.loc[multi_domain['FA-DOMID']==domain2, 'FA-UNIID'].item()
    
    dict_2domains['d1_start'] = int(multi_domain.loc[multi_domain['FA-DOMID']==domain1, 'FA-PDBREG-START'].item())
    dict_2domains['d1_start_unp'] = int(multi_domain.loc[multi_domain['FA-DOMID']==domain1, 'FA-UNIREG-START'].item())
    
    # d1_end = multi_domain.loc[multi_domain['FA-DOMID'], 'FA-PDBREG-END']
    dict_2domains['d1_end'] = loop['start'] - 1
    dict_2domains['d1_end_unp'] = loop['start_unp'] - 1
    
    # d2_start = multi_domain.loc[multi_domain['FA-DOMID'], 'FA-PDBREG-START']
    dict_2domains['d2_start'] = loop['end'] + 1
    dict_2domains['d2_start_unp'] = loop['end_unp'] + 1
    dict_2domains['d2_end'] = int(multi_domain.loc[multi_domain['FA-DOMID']==domain2, 'FA-PDBREG-END'].item())
    dict_2domains['d2_end_unp'] = int(multi_domain.loc[multi_domain['FA-DOMID']==domain2, 'FA-UNIREG-END'].item())
    
    return dict_2domains
    
def save_missingPercentage(dict_domain_loops, path_dssp_full_chainid, path_loops_domain_smoothSS_2_5_perMiss):
    '''
    Calculate the percentage of missing residues for each loop.
    params:
        dict_domain_loops - dict, domain loops
        path_dssp_full_chainid - chain (sequence) information including missing residues
        path_loops_domain_smoothSS_2_5_perMiss - save to file path
    '''
    for pid_chain, values in dict_domain_loops.items():
        print(pid_chain)
        # load chain table
        df_chain = pd.read_json(os.path.join(path_dssp_full_chainid, pid_chain+'.json'))
        for i in range(len(values)):
            loop = values[i]
            start = loop['start']
            end = loop['end']
            length = loop['length']
            
            # convert auth_seqidx from float into int. 
            # df_chain['auth_seqidx'] = [int(x) for x in df_chain['auth_seqidx']]
            
            # get loop region and calculate missing percentage
            df_loop = df_chain[df_chain['auth_seqidx'].isin(range(start, end+1))]
            len_missing = sum(df_loop['missing'])
            per_missing = len_missing/length

            # save missing percentage to dictionary
            dict_domain_loops[pid_chain][i]['miss_length'] = len_missing
            dict_domain_loops[pid_chain][i]['miss_percentage'] = per_missing
    dump_dict2json(dict_domain_loops, path_loops_domain_smoothSS_2_5_perMiss)