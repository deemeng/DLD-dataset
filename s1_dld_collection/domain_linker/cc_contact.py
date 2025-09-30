import os
import pandas as pd
import numpy as np

from Bio.PDB import MMCIFParser, FastMMCIFParser
from utils.common import load_tab, save_tab, dump_dicts2jsons, lower_pdbid, read_json2dict

'''
Functions related to calcute distance between 2 CAs from 2 adjacent Domains.

Di
Mar, 19, 2023

'''

def getCoords(pdbid: str, chainid: str, path_cif) -> dict:
    '''
    get coordinates of all residues from a pdbid.
    
    params:
        pdbid - PDBID
        path_cif - mmcif file path (.cif)
    
    return:
        seqidx_coords - dict, {auth_seq_idx, [x, y, z]}
    '''
    seqidx_coords = {}

    cif_parser = FastMMCIFParser()
    structure = cif_parser.get_structure(pdbid, path_cif)

    model = structure[0] 
    chain = model[chainid.upper()]
    for residue in chain.get_residues():
        try:
            seq_idx = residue.get_id()[1]
            seqidx_coords[seq_idx] = residue['CA'].coord
        except:
            break
    return seqidx_coords

def cal_distance(coord1, coord2) -> float:
    '''
    Calcaulate Euclidean distance between 2 points.
    '''
    # the Euclidean distance is the l2 norm, and the default value of the ord parameter in numpy.linalg.norm is 2
    distance = np.linalg.norm(coord1 - coord2)
    return distance

def contact(dict_2domains: dict, seqidx_coords:dict, t=5, log=None):
    '''
    we assumed an interdomain CA-CA contact when the distance between the side-chain carbon atoms of three or more hydrophobic residues was less than 4.663 A ̊ 
    params:
        dict_2domains - {'d1_start', 'd1_end', 'd2_start', 'd2_end'}
        seqidx_coords - dict, {auth_seq_idx, [x, y, z]}
        t - the threshold of distance between 2 CA, 5AS
        
    return:
        residue_distance[f'{r1}_{r2}']
    '''
    
    residue_distance = {}
    
    for r1 in range(dict_2domains['d1_start'], dict_2domains['d1_end']+1):
        for r2 in range(dict_2domains['d2_start'], dict_2domains['d2_end']+1):
            if (r1 in seqidx_coords) & ((r2 in seqidx_coords)):
                dist = cal_distance(seqidx_coords[r1], seqidx_coords[r2])
                if dist < t:
                    # print(f'\ndistance: {dist}\nresidues: \nr1: {r1} - {seqidx_coords[r1]}, r2: {r2} - {seqidx_coords[r2]}')
                    residue_distance[f'{r1}_{r2}'] = str(dist)
                    if log is not None:
                        log.info(f'\ndistance: {dist}\nresidues: \nr1: {r1} - {seqidx_coords[r1]}, r2: {r2} - {seqidx_coords[r2]}')
                        
    return residue_distance