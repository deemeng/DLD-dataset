from Bio.PDB import MMCIFParser, is_aa
import numpy as np
from params import *

from utils.common import load_tab, read_json2dict
from domain_linker.extract_linkers import get_2domain_start_end


'''
This sript aims to rpovide functions to calculate the side-chain carbons distance between two adjicent domains. 
We only record the distance less than 5A. 

Functions & Globular parameters
'''

# Parameters
CUTOFF = 5  # Ångström cutoff for side-chain contact

# List of hydrophobic residues
hydrophobic_residues = ['VAL', 'LEU', 'ILE', 'PHE', 'PRO', 'MET', 'TRP', 'TYR']

# Helper: is side-chain carbon?
def is_sidechain_carbon(atom):
    return atom.element == 'C' and atom.get_id() not in ('C', 'CA')

# Collect side-chain carbons for hydrophobic residues
def collect_sidechain_carbons(path_cif, chainid, domain_range):
    # Parse structure
    parser = MMCIFParser()
    structure = parser.get_structure("struct", path_cif)
    model = structure[0]
    chain = model[chainid.upper()]
    atoms = []
    
    for residue in chain:
        if not is_aa(residue):
            continue
        res_name = residue.get_resname()
        res_id = residue.get_id()[1]
        if res_id in domain_range and res_name in hydrophobic_residues:
            for atom in residue:
                if is_sidechain_carbon(atom):
                    atoms.append((res_id, atom))
    return atoms

def calcul_dist(domain1_atoms, domain2_atoms):
    # Compute interdomain contacts
    contacts = []
    for res1, atom1 in domain1_atoms:
        coord1 = atom1.get_coord()
        for res2, atom2 in domain2_atoms:
            coord2 = atom2.get_coord()
            dist = np.linalg.norm(coord1 - coord2)
            if dist < CUTOFF:
                contacts.append((res1, res2, atom1.get_id(), atom2.get_id(), str(round(dist, 2))))
    
    return contacts