import os
import pandas as pd

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBList

from utils.common import load_tab
from loop.pdb_files import get_PDBfiles


ROOT = os.path.realpath('..')

path_recollect = os.path.join(ROOT, 'domain_linker/recollect_dsAll/')

# folder
path_pdb_cif = os.path.join(path_recollect, 'pdb/mmcif')

# multi_domain file path
path_scop_folder = os.path.join(path_recollect, 'scop2')
path_tab_scop_FA_continuous_uni_multi_domian = os.path.join(path_scop_folder, 'scop-cla-latest-table-FAcontinuousUniMultiDomain.tsv')

# download pdbfiles
list_pdbid = load_tab(path_tab_scop_FA_continuous_uni_multi_domian)['FA-PDBID'].unique()
get_PDBfiles(list_pdbid, path_pdb_cif)