import os
import pandas as pd

def smothing_ss(sec_structure_3: list, coil=True, threshold=4):
    '''
    coil=False - Set SS(H, E) to C if the len(SS)<threshold (4)
    coil=True - Set C to H if the len(C)<threshold (4)
    
    params:
        sec_structure_3 - a list of secondary structures, H, E, C
        coil - True  => smoothing coil
               False => smoothing H, E
        threshold - the minimum length of a loop or a SS
    return:
        Smoothing sec_structure_3
        corresponding marks
    '''
    if coil:
        ss_mark = [1 if ss=='C' else 0 for ss in sec_structure_3]
    else:
        # C: 0, E:2, H:1
        ss_mark = [0 if ss=='C' else 2 if ss=='E' else 1 for ss in sec_structure_3]
    
    count = 0
    for i in range(len(ss_mark)):
        if ss_mark[i] == 0:
            if count < threshold:
                ss_mark[i-count: i] = [0] * count
                if coil:
                    sec_structure_3[i-count: i] = ['H'] * count
                else:
                    if 'E' not in sec_structure_3[i-count: i]:
                        sec_structure_3[i-count: i] = ['C'] * count
            count = 0
        else:
            count += 1

    # ✅ handle leftover at the end (C-terminal)
    if count > 0 and count < threshold:
        ss_mark[len(ss_mark)-count:] = [0] * count
        if coil:
            sec_structure_3[len(ss_mark)-count:] = ['H'] * count
        else:
            if 'E' not in sec_structure_3[len(ss_mark)-count:]:
                sec_structure_3[len(ss_mark)-count:] = ['C'] * count

    return sec_structure_3, ss_mark

# currently don't use this function anymore.
# missing/unmodeled residues are disordered.
def check_unmodeled(loop: dict):
    '''
    Check if a region includes un modeled residues or not;
    len(seq) == end - start + 1
    
    params:
        loop - {'start', 'end', 'seq'}
    return:
        True - include unmodeled residue
        False - no unmodeled residue
    '''
    return len(loop['seq']) < (int(loop['end']) - int(loop['start']) + 1)

def get_2steps_smoothing_ss(df_chain: pd.DataFrame, t1=3, t2=4):
    '''
    Smoothing the SS elements of sec_structure_3 by length 2 (len_SS<3, set SS to coil), save the result to sec_structure_3_s2 column.
    Then, smoothing the SS elements of sec_structure_3_2 by length 3 (len_SS<4, set SS to coil), save the result to sec_structure_3_s2s3 column.
    
    params:
        df_chain - current chain information
        t1 - int, the threshold for first step of smoothing
        t2 - int, the threshold for second step of smoothing
    return:
        df_chain - chain information with 2 more new columns, sec_structure_3_s2 and sec_structure_3_s2s3
    '''
    # column names
    sec_structure_3_first = 'sec_structure_3_s'+str(t1-1)
    sec_structure_3_second = 'sec_structure_3_s'+str(t1-1)+'s'+str(t2-1)
    df_chain[sec_structure_3_first], _ = smothing_ss(list(df_chain['sec_structure_3']), coil=False, threshold=t1)
    df_chain[sec_structure_3_second], _ = smothing_ss(list(df_chain[sec_structure_3_first]), coil=False, threshold=t2)
    
    return df_chain

def get_loops4one_chain(df_chain: pd.DataFrame, ss_1=False, ss_2=False, t1=3, t2=4) -> list[dict]:
    '''
    get all loops in a chain.
    
    params:
        df_chain - current chain information
        ss_1 - False: do not use sec_structure_3_s2
               True: use sec_structure_3_s2
        ss_2 - False: do not us sec_structure_3_s2s3
               True: use sec_structure_3_s2s3
        t1 - int, the threshold for first step of smoothing
        t2 - int, the threshold for second step of smoothing
    return:
        loops4one_chain - [{'start': , 'end': , 'loop_seq': ,}, ]
    '''
    sec_structure_3_first = 'sec_structure_3_s'+str(t1-1)
    sec_structure_3_second = 'sec_structure_3_s'+str(t1-1)+'s'+str(t2-1)
    
    if ss_1:
        df_chain = get_2steps_smoothing_ss(df_chain)
        col_name = sec_structure_3_first
    elif ss_2:
        df_chain = get_2steps_smoothing_ss(df_chain)
        col_name = sec_structure_3_second
    else:
        col_name = 'sec_structure_3'
    
    '''
    columns in df_chain:
    
    ['dssp_idx', 'pdbid', 'sequence', 'chain', 'auth_seqidx',
       'sec_structure_8', 'sec_structure_3', 'NH2O_1_relidx', 'NH2O_1_energy',
       'O2NH_1_relidx', 'O2NH_1_energy', 'NH2O_2_relidx', 'NH2O_2_energy',
       'O2NH_2_relidx', 'O2NH_2_energy', 'missing', 'asym_id', 'seq_id',
       'unp_res', 'unp_num', 'unp_acc', 'unp_segment_id', 'unp_instance_id']
    '''
    sec_structure_3 = list(df_chain[col_name])
    sequence = list(df_chain['sequence'])
    pdb_seqidx = list(df_chain['auth_seqidx'])
    missing = list(df_chain['missing'])
    dssp_key_str = list(df_chain['dssp_key_str'])
    ####
    # uniprot information
    ####
    
    unp_res = list(df_chain['unp_res']) # residue
    unp_num = list(df_chain['unp_num']) # index
    unp_acc = list(df_chain['unp_acc']) # accession number
    
    # smoothing the SS
    ss_smooth, ss_mark = smothing_ss(sec_structure_3)
    
    loops4one_chain = []
    count = 0
    for i in range(len(ss_mark)):
        if ss_mark[i] == 0 or i==len(ss_mark)-1:
            if count != 0:
                if i==len(ss_mark)-1: # count last loop
                    i = i+1
                    count = count + 1
                    
                loop = {}
                # auth_idx
                loop['start'] = pdb_seqidx[i - count]
                loop['end'] = pdb_seqidx[i - 1]
                loop['seq_id'] = pdb_seqidx[i-count: i]
                
                loop['start_unp'] = unp_num[i - count]
                loop['end_unp'] = unp_num[i - 1]
                loop['seq_id_unp'] = unp_num[i-count: i]
                
                # auth_seq
                loop['seq'] = sequence[i-count: i]
                loop['seq_unp'] = unp_res[i-count: i]
                loop['dssp_key_str'] = dssp_key_str[i-count: i]
                
                # missing
                loop['missing'] = missing[i-count: i]
                loop['miss_length'] = sum(missing[i-count: i])
                loop['miss_percentage'] = sum(missing[i-count: i])/count
                
                loop['unp_acc'] = unp_acc[0]
                loop['length'] = count
                
                '''
                # if do not include unmodeled residues
                if not check_unmodeled(loop):
                    loops4one_chain.append(loop)
                '''
                # unmodeled/missing residues are disordered
                loops4one_chain.append(loop)
                # reset count
                count = 0
            else:
                continue
        else:
            count = count + 1
            
    return loops4one_chain

def get_loops(list_pidChain: list, path_dssp_chainid: str, ss_1=False, ss_2=False, t1=3, t2=4):
    '''
    Get all loops in All chains
    
    params:
        list_pidChain - [pdbid_chainids]
        path_dssp_chainid - chain json files folder
        ss_1 - False: do not us sec_structure_3_s2
               True: use sec_structure_3_s2
        ss_2 - False: do not us sec_structure_3_s2s3
               True: use sec_structure_3_s2s3
        t1 - int, the threshold for first step of smoothing
        t2 - int, the threshold for second step of smoothing
    return:
        loops - {pdbid_chainid: [loops in the current chain]}
    '''
    loops = {}
    for pid_chain in list_pidChain:
        # load JSON file to DataFrame
        # print(os.path.join(path_dssp_chainid, pid_chain+'.json'))
        df_chain = pd.read_json(os.path.join(path_dssp_chainid, pid_chain+'.json'))
        loops[pid_chain] = get_loops4one_chain(df_chain, ss_1, ss_2, t1, t2)
        
    return loops

def check_terminus (loop: dict, domains_in_chain: list):
    '''
    check if the loop is a terminu based on all domain regions in a chain.
    
    params:
        loop - {start:, end:, seq:, len:,}
        domains_in_chain - [[domain1_start, domain1_end], [domain2_start, domain2_end], ...]
    return:
        dict, {'terminu': , 'domain':}
            is-terminu - True, terminu
                      False, intra-domain loops
            domain - DomainID
                  
    '''
    
    terminu = True
    
    l_start = loop['start']
    l_end = loop['end']
    
    domain_count = 0
    
    for domain in domains_in_chain:
        d_start = domain[0]
        d_end = domain[1]
        d_id = domain[2]
        
        if (l_start<d_start) and (domain_count==0):
            # print(domain, l_start, l_end, 'N')
            # N-terminu
            return {'terminu': terminu, 'domain':d_id, 'tag': 'N'}
            
        if (l_start>=d_start) and (l_end<=d_end):
            # print(domain, l_start, l_end)
            terminu = False
            # intra-domain loop
            return {'terminu': terminu, 'domain':d_id, 'tag': 'intra'}
        
        if (l_end>=d_end) and (domain_count==(len(domains_in_chain)-1)):
            # print(domain, l_start, l_end, 'C')
            return {'terminu': terminu, 'domain':d_id, 'tag': 'C'}
        
        # loops overlap with one domain: the end part of the domain
        if ((l_start>=d_start) & (l_start<=d_end)):
            # print(domain, l_start, l_end, 'overlap1_end')
            return {'terminu': False, 'domain':d_id, 'tag': 'overlap1_end'}
        
        # loops overlap with one domain: the start part of the domain
        if ((l_end>=d_start) & (l_end<=d_end)):
            # print(domain, l_start, l_end, 'overlap1_start')
            return {'terminu': False, 'domain':d_id, 'tag': 'overlap1_start'}
        
        domain_count += 1
        
    # print([], l_start, l_end, 'overlap0')
    
    # loops do not overlap with any domain
    return {'terminu': False, 'domain':-1, 'tag': 'overlap0'}


def check_overlap_2domain (loop: dict, domains_in_chain: list):
    '''
    check if a loop is overlap with 2 domains. 1step smoothing loop. For irregular examples, 3zx8_a.
    [check if this 'Intra-domain-loop' is overlapping with 2 domains.]
    After 2step of smoothing, one of the doamin becomes a loop.
    
    params:
        loop - {start:, end:, seq:, len:,}
        domains_in_chain - [[domain1_start, domain1_end, domainID], [domain2_start, domain2_end, domainID], ...]
    return:
        overlap_domain - {overlap, domains: []}
            overlap: True, overlap with two domains
                      False, do not overlap with 2 domains
            domains: The two domain it overlapped with.
    '''
    ### concat all starts & ends. then check it.!!!! think about the terminus
    overlap = False
    
    l_start = loop['start']
    l_end = loop['end']

    for d in range(len(domains_in_chain)-1):
        d1_start = domains_in_chain[d][0]
        d1_end = domains_in_chain[d][1]
        d2_start = domains_in_chain[d+1][0]
        d2_end = domains_in_chain[d+1][1]
        if (l_start in range(d1_start,d1_end+1)) & (l_end in range(d2_start, d2_end+1)):
            overlap = True
            return {'overlap': overlap, 'domains': [domains_in_chain[d][2], domains_in_chain[d+1][2]]}
    return {'overlap': overlap}