import numpy as np

'''
show loops_info related functions.

Di
Mar, 20, 2023
'''


def count_loops(dict_loops):
    '''
    count numbers of loops, chains, and pdb-entries
    
    params:
        dict_loops - {pdbid_chainid: [loops]}
    return:
        count - num of loops
        count_chain - num of chains
        count_pdbid - num of PDB entries
    '''
    count = 0
    count_chain = 0
    count_pdbid = 0
    
    list_pdbid = []
    
    for key, loop in dict_loops.items():
        count = count + len(loop)
        count_chain = count_chain + 1
        list_pdbid.append(key.split('_')[0])
    count_pdbid = len(set(list_pdbid))
    
    return count, count_chain, count_pdbid

def print_loop_count(dict_loops):
    '''
    print loop count info
    '''
    count, count_chain, count_pdbid = count_loops(dict_loops)
    print(f'The  number of Loops: {count}')
    print(f'The  number of Chains: {count_chain}')
    print(f'The  number of PDB entries: {count_pdbid}')