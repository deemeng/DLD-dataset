import os

ROOT = os.path.realpath('..')

path_recollect = os.path.join(ROOT, 'data/recollect_dsAll/')
path_orgDSAll = os.path.join(ROOT, 'data/origin_dsAll/')

path_dsAll = os.path.join(path_orgDSAll, 'DS-All_LinkerList.txt')
path_dsAll_tab = os.path.join(path_orgDSAll, 'DS-All.tab')
# folder
path_pdb_cif = os.path.join(path_recollect, 'pdb/mmcif')
path_pdbe_cif = os.path.join(path_recollect, 'pdbe/updated_mmcif')
path_pdbe_uniprot_pdbid = os.path.join(path_recollect, 'pdbe/uniprot_pdbid')
path_dssp_pdbid = os.path.join(path_recollect, 'dssp/pdbid')
path_dssp_chainid = os.path.join(path_recollect, 'dssp/chainid')
path_dssp_full_pdbid = os.path.join(path_recollect, 'dssp/full_pdbid')
# path_dssp_full_uniprot_pdbid = os.path.join(path_recollect, 'dssp/full_uniprot_pdbid')
path_dssp_full_uniprot_chainid = os.path.join(path_recollect, 'dssp/full_uniprot_chainid')

path_loops = os.path.join(path_recollect, 'loops')

# files
path_loops_all = os.path.join(path_loops, 'all/loops.json')
path_loops_all_smoothSS_first = os.path.join(path_loops, 'all/loops_smoothSS_first.json')
path_loops_all_smoothSS_second = os.path.join(path_loops, 'all/loops_smoothSS_second.json')
path_loops_domain_smoothSS_second = os.path.join(path_loops, 'all/domain_loops_smoothSS_second.json')

path_loops_domain_perMiss = os.path.join(path_loops, 'all/domain_loops_perMiss.json')
path_loops_domain_smoothSS_first_perMiss = os.path.join(path_loops, 'all/domain_loops_smoothSS_first_perMiss.json')
path_loops_domain_smoothSS_second_perMiss = os.path.join(path_loops, 'all/domain_loops_smoothSS_second_perMiss.json')

# excepltions
path_dssp_error = os.path.join(path_recollect, 'dssp/parse_error.txt')
path_nomissing_residue = os.path.join(path_recollect, 'dssp/no_missing_residue_pdbids.txt')

# multi_domain file path
path_scop_folder = os.path.join(path_recollect, 'scop2')
path_tab_scop_FA_continuous_uni_multi_domian = os.path.join(path_scop_folder, 'scop-cla-latest-table-FAcontinuousUniMultiDomain.tsv')
path_tab_scop_FA_continuous_uni_multi_domian_withResolution = os.path.join(path_scop_folder, 'scop-cla-latest-table-FAcontinuousUniMultiDomainWithResolution .tsv')
path_tab_scop_FA_continuous_uni_multi_domian_lenResolution_filtered = os.path.join(path_scop_folder, 'scop-cla-latest-table-FAcontinuousUniMultiDomain_lenResolution_filtered.tsv')
path_tab_scop_FA_continuous_uni_multi_parse_domian = os.path.join(path_scop_folder, 'scop-cla-latest-table-FAcontinuousUniMultiParseDomain.tsv')
path_tab_scop_FA_continuous_uni_multi_parse_domian_dist = os.path.join(path_scop_folder, 'scop-cla-latest-table-FAcontinuousUniMultiParseDomainDist.tsv')

# domain_dataset, without sequence infomation
path_domain_dataset_pre = os.path.join(path_loops, 'all/domain_dataset_pre.json')

# final

# N-/C-terminus
path_terminus = os.path.join(path_loops, 'final/terminus.json')
path_terminus_tab = os.path.join(path_loops, 'final/terminus.tsv')

# pdb_chain_unp
path_pdb_chain_unp = os.path.join(path_loops, 'final/pdb_chain_unp.json')
path_idxMap_pdb_chain_unp = os.path.join(path_loops, 'final/mapIdx_pdb_chain_unp.tsv')

# inter-domain-linkers/Loops
path_loops_domain_dist_type = os.path.join(path_loops, 'all/domain_loops_dist_type.tsv')

path_loops_domain_hbond_contact_smoothSS_second = os.path.join(path_loops, 'final/all_domain_linker_hbond_contact_smoothSS_second.json')

path_loops_domain_dist_hbond_contact_smoothSS_second = os.path.join(path_loops, 'final/all_domain_linker_dist_hbond_contact_smoothSS_second.tsv')

path_independent_domain_linker = os.path.join(path_loops, 'final/independent_domain_linker.tsv')
path_dependent_domain_linker = os.path.join(path_loops, 'final/dependent_domain_linker.tsv')

path_independent_domain_linker_json = os.path.join(path_loops, 'final/independent_domain_linker.json')
path_dependent_domain_linker_json = os.path.join(path_loops, 'final/dependent_domain_linker.json')
# intra-domain
# intra-domain loops (all loops = intra-domain_smooSS_2 loops + inter-domain_smooSS_2_3 loops)
path_loops_intra_domain_smoothSS_first_2step = os.path.join(path_loops, 'final/intra_domain_loops_smoothSS_first_2step.json')
path_loops_intra_domain_smoothSS_first_2step_tab = os.path.join(path_loops, 'final/intra_domain_loops_smoothSS_first_2step.tsv')

# domain_dataset
path_domain_dataset = os.path.join(path_loops, 'final/domain_dataset.json')
path_domain_dataset_tab = os.path.join(path_loops, 'final/domain_dataset.tsv')