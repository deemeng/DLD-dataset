import os
import pandas as pd

from utils.common import save_tab, load_tab

from loop.extract_loop import get_loops
from params import *

'''
Save parsable pdb entries to a new table

15 Mar, 2023
Di
'''

# load files
multi_domain_all = load_tab(path_tab_scop_FA_continuous_uni_multi_domian_lenResolution_filtered)
list_pdbid = multi_domain_all['FA-PDBID'].unique()
list_pdbid_error = pd.read_csv(path_dssp_error, names=['chainid', 'error'])['chainid'].unique()

# keep PDB entries without errors
list_pdbid_error = [x.upper() for x in list_pdbid_error]
pdbids = list(set(list_pdbid) - set(list_pdbid_error))
multi_domain = multi_domain_all[multi_domain_all['FA-PDBID'].isin(pdbids)]

# delete domains with unbigious start & end, e.g. 1p 100p
multi_domain = multi_domain[[True if x.isnumeric() else False for x in multi_domain['FA-PDBREG-START']]]
multi_domain = multi_domain[[True if x.isnumeric() else False for x in multi_domain['FA-PDBREG-END']]]

# delete chains only have one domain available
multi_domain = multi_domain.groupby(['FA-PDBID', 'FA-CHAINID']).filter(lambda x: len(x)>1).sort_values(['FA-PDBID', 'FA-CHAINID', 'FA-PDBREG-START'])

# save
save_tab(multi_domain, path_tab_scop_FA_continuous_uni_multi_parse_domian)