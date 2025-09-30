import os
import pandas as pd
import numpy as np

from Bio.PDB import MMCIFParser, FastMMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.DSSP import DSSP

from loop.pdb_files import get_PDBfiles, parse_mmcif
from utils.common import load_tab, save_tab, dump_dicts2jsons, dump_dict2json, read_json2dict
from utils.common import int_domainREG

from domain_linker.cc_contact import contact, getCoords
from domain_linker.extract_linkers import get_2domain_start_end
from domain_linker.hbond import check_2Hbond
from params import *

from utils.log import logger
'''
Aim to separate domain linkers to
    - inter-domain linkers
    - intra-domain linkers: 1. Hydrophobic interaction: 3 or more dist(CA-from-two-domains)<4.663
                            2. three or more H-bond between two domains
    !!!! check Hbonds using path_dssp_chainid, not path_dssp_full_chainid. Because path_dssp_full_chainid includes missing residues
    !!!! Domain loops were extracted from smoothSS_2_3

May, 3, 2023.
Di
'''

# log message
log = logger('cc_contact_smoothSS_second.log')

# loadin domain info
multi_domain = load_tab(path_tab_scop_FA_continuous_uni_multi_parse_domian)

# convert non-numeric values to numeric values
# col_names = ['FA-PDBREG-START', 'FA-PDBREG-END']
# multi_domain = int_domainREG(multi_domain, col_names)

# load in all domain-loops, or we can call it domain-linkers including both inter-domain & intra-domain linkers
try: 
    print('path_loops_domain_perMiss is available!')
    log.info('path_loops_domain_perMiss is available!')
    dict_domain_loops = read_json2dict(path_loops_domain_smoothSS_second_perMiss)
except:
    print('path_loops_domain!')
    log.info('path_loops_domain!')
    dict_domain_loops = read_json2dict(path_loops_domain_smoothSS_second)

# as the domain-linker index
loop_count = 0

# save all domain linkers
domain_linkers = {}

# indicate domian-linkers are inter or intra
# {'inter': [], 'intra': []}
inter_intra_linkers = {}

###
# Loop through all domain-linkers and separate them 
###

for pid_chain_unp, domain_loop in dict_domain_loops.items():
    print(pid_chain_unp)
    
    pdbid = pid_chain_unp.split('_')[0]
    chainid = pid_chain_unp.split('_')[1]
    unp = pid_chain_unp.split('_')[2]
    # domain_loop = dict_domain_loops[pid_chain_unp]
    path_cif = os.path.join(path_pdb_cif, pdbid+'.cif')
    
    for loop in domain_loop:
        # loop['loop_id'] = loop_count
        loop['pdbid'] = pdbid
        loop['chainid'] = chainid
        
        # 1. check hydrophobic interactions
        # get coordinations of residues
        seqidx_coords = getCoords(pdbid, chainid, path_cif)
        dict_2domains = get_2domain_start_end(loop, multi_domain)
        residue_distance = contact(dict_2domains, seqidx_coords, t=5, log=log)
        
        # 2. H-bond
        # check if there are more than 2 Hbond between 2 domains
        df_chain = pd.read_json(os.path.join(path_dssp_chainid, f'{pdbid.lower()}_{chainid.lower()}.json'))
        hbond_idx = check_2Hbond(df_chain, loop, dict_2domains)
        
        loop['hbonds'] = hbond_idx
        loop['c_c_contacts'] = residue_distance
        
        domain_linkers[loop_count] = loop
        loop_count = loop_count + 1
###
# Save domain linkers and the separate (inter or intra) information
###
print('Saving Linkers ... ')

dump_dict2json(domain_linkers, path_loops_domain_hbond_contact_smoothSS_second)
print('Done!')