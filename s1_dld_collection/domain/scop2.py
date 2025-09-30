import os
import pandas as pd
import numpy as np


###
# SCOP: Domain data collection
# version: SCOP2
# date: 11 Mar, 2023
# Functions related to SCOP data processing
###

# 1. load scop file to pd.Dataframe
def scop2df(path_scop: str) -> pd.DataFrame:
    '''
    params:
        path_scop - file path of scop cla file.
    return:
        df_scop - load scop cla file to pandas Dataframe.
    '''
    
    # a. load column names

    # open file & load in lines
    file = open(path_scop)
    content = file.readlines()
    # read the 6th line
    column_names = content[5][2:].split()

    # b. load data to dataframe
    
    df_scop = pd.read_csv(path_scop, skiprows=6, names=column_names, sep=' ')
    print(f'length of df_scop: {df_scop.shape[0]}')
    
    return df_scop

# 2. Separate SCOPCLA column
def scopcla(df_scop: pd.DataFrame) -> pd.DataFrame:
    '''
    params:
        df_scop - load from scop file
    return:
        df_scopcla - A dataframe with columns extra 'TP', 'CL', 'CF', 'SF', 'FA'.
    '''
    df_cla = pd.DataFrame([list(map(lambda x: x[3:], x.split(','))) for x in df_scop['SCOPCLA']], columns=['TP', 'CL', 'CF', 'SF', 'FA'])
    df_scopcla = pd.concat([df_scop, df_cla], axis=1)
    return df_scopcla


# 3. split region field into individual start-end
def separate_UNIreg(df_scop, level='FA'):
    '''
    params:
        df_scop - Dataframe datatype, including domain information from SCOP.
        level - FA or SF
        
    return:
        A dataframe with DOMID and the starts and ends of all region/domain.
    '''
    list_domains = []
    
    for domainID, reg in zip(df_scop[level+'-DOMID'], df_scop[level+'-UNIREG']):
        for r in reg.split(','):
            se = r.split('-')
            s = int(se[0])
            e = int(se[1])

            # FA domains
            dict_domain = {}
            dict_domain[level+'-DOMID'] = domainID

            dict_domain[level+'-UNIREG-START'] = s
            dict_domain[level+'-UNIREG-END'] = e

            list_domains.append(dict_domain)
    return pd.DataFrame(list_domains)

# 4. split region field into individual start-end
def separate_PDBreg(df_scop, level='FA'):
    '''
    params:
        df_scop - Dataframe datatype, including domain information from SCOP.
        level - FA or SF
        
    return:
        A dataframe with DOMID and the starts and ends of all region/domain.
    Note: start-end are not integers since the following situation exist.
        FA-DOMID	FA-PDBID	FA-PDBREG	FA-UNIID	FA-UNIREG	FA: 
        8054448	1KJ2	B:1-116A	P04214	22-136	4007530
    '''
    list_domains = []
    
    for domainID, reg in zip(df_scop[level+'-DOMID'], df_scop[level+'-PDBREG']):
        for r in reg.split(','):
            list_r = r.split(':')
            chain = list_r[0]
            se = list_r[1].split('-')
            try: 
                if len(se)==2:
                    # FA-DOMID	FA-PDBID	FA-PDBREG	FA-UNIID	FA-UNIREG	FA: 
                    # 8054448	1KJ2	B:1-116A	P04214	22-136	4007530
                    s = se[0]
                    e = se[1]
                else:
                    # situation
                    # FA-DOMID	FA-PDBID	FA-PDBREG	FA-UNIID	FA-UNIREG	FA: 
                    # 8094978	5WQ0	D:-3-125	A0A1W2VN01	1-129	4003632
                    s = '-' + se[1]
                    e = se[2]
            except:
                print(f'! The function does not work for: {domainID}')

            # FA domains
            dict_domain = {}
            dict_domain[level+'-DOMID'] = domainID

            dict_domain[level+'-PDBREG-START'] = s
            dict_domain[level+'-PDBREG-END'] = e
            dict_domain[level+'-CHAINID'] = chain

            list_domains.append(dict_domain)
    return pd.DataFrame(list_domains)

def get_domain_start_end(multi_domain: pd.DataFrame) -> dict:
    '''
    Get all domains' start and end in every PDB chain. Return {(pdbid, chainid): [[domain1_start, domain1_end], []]}
    
    params:
        multi_domain - domain information from scop
    return:
        dict_chain_domainID - {(pdbid, chainid): [[domain1_start, domain1_end, domian1_ID], [domain2_start, domain2_end, domain2_ID], ...], }
    '''
    multi_domain = multi_domain.sort_values(by=['FA-PDBID', 'FA-CHAINID', 'FA-PDBREG-START'])
    dict_chain_domainID = multi_domain.groupby(['FA-PDBID', 'FA-CHAINID'])[['FA-PDBREG-START', 'FA-PDBREG-END', 'FA-DOMID']].apply(lambda g: g.values.tolist()).to_dict()
    return dict_chain_domainID

def get_domain_dataset(multi_domain: pd.DataFrame) -> dict:
    '''
    Get the original domain dataset from scop, save it to a dict.
    
    params:
        multi_domain - domain information from scop
    return:
        dict_domain_dataset - {'pdb_chainid': , 'start': , 'end': , 'start_unp': , 'end_unp': , 'unp_acc': ,}
    
    '''
    dict_domain_dataset = {}
    
    for idx, row in multi_domain.iterrows():
        dict_domain_dataset[row['FA-DOMID']] = {'pdb_chainid': row['FA-PDBID'].lower()+'_'+row['FA-CHAINID'].lower(), 
                                                      'start': row['FA-PDBREG-START'], 'end': row['FA-PDBREG-END'],
                                                     'start_unp': row['FA-UNIREG-START'], 'end_unp': row['FA-UNIREG-END'],
                                                     'unp_acc': row['FA-UNIID'],
                                                     }
    return dict_domain_dataset