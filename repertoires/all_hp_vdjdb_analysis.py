#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 13:49:43 2024

@author: wezi
"""

import pandas as pd
import scipy.sparse
from tcrdist.breadth import get_safe_chunk
from tcrdist.repertoire import TCRrep
from tcrdist.join import join_by_dist
from tcrdist.background import sample_britanova
from tcrdist.background import get_stratified_gene_usage_frequency
from tcrdist.pgen import OlgaModel

CPUS = 4

# load up the repertoire
file = 'all_hp_bal_tcrb.tsv'

df_bulk = pd.read_csv('all_hp_bal_tcrb.tsv', sep='\t')

# remove the rows where the CDR3 is incompletely sequenced
filter = df_bulk['cdr3_b_aa'].str.contains('_')
df_bulk = df_bulk[~filter]

# sort by count (allows distinction between identical aa with differing nt
df_bulk = df_bulk.sort_values('count', ascending = False).reset_index(drop = True)

df_bulk['rank'] = df_bulk.index.to_list()

tr_bulk = TCRrep(
    cell_df = df_bulk, 
    organism = 'human', 
    chains = ['beta'], 
    db_file = 'alphabeta_gammadelta_db.tsv', 
    compute_distances = True)

tr_bulk.compute_distances()

# Weird bug makes clone ID one lower when referred to in nn list. Fix
tr_bulk.clone_df['clone_id'] = tr_bulk.clone_df['clone_id'] - 1

# Sixfold weighting for CDR3 in TCRdist units
tr_bulk.weights_b = {'cdr3_b_aa': 6, 'pmhc_b_aa': 1, 'cdr2_b_aa': 1, 'cdr1_b_aa': 1}
chunk_size = get_safe_chunk(tr_bulk.clone_df.shape[0], tr_bulk.clone_df.shape[0])
tr_bulk.compute_sparse_rect_distances(radius = 72, chunk_size = chunk_size)
scipy.sparse.save_npz(f'{file}.tr_bulk.rw_beta_csrmat.npz', tr_bulk.rw_beta)

tr_bulk.rw_beta=scipy.sparse.load_npz(f'{file}.tr_bulk.rw_beta_csrmat.npz')

# load up cord blood sample. Easier to do this even if not needed in this
# analysis
df_cord = sample_britanova(960000, random_state=1)

ts = get_stratified_gene_usage_frequency(replace = True)

# Uncomment above and change cell_df back to df_nn if you want to filter by
# clones with > k_nn neighbours
tr_nn = TCRrep(
    cell_df = tr_bulk.clone_df,
    organism = "human",
    chains = ['beta'],
    deduplicate = False,
    compute_distances = False)

tr_nn.cpus = CPUS

# Sixfold weighting for CDR3 in TCRdist units
tr_nn.weights_b = {'cdr3_b_aa': 6, 'pmhc_b_aa': 1, 'cdr2_b_aa': 1, 'cdr1_b_aa': 1}

tr_cord = TCRrep(
    cell_df = df_cord,
    organism = "human",
    chains = ['beta'],
    cpus = CPUS,
    compute_distances = False)

# Have a look to see if clones have higher within-repertoire neighbours
# compared to neighbours within cord blood (which would show antigenic
# selection)
chunk_size = get_safe_chunk(tr_nn.clone_df.shape[0], tr_cord.clone_df.shape[0])

tr_nn.compute_sparse_rect_distances(df = tr_nn.clone_df, df2 = tr_cord.clone_df, radius = 72, chunk_size=chunk_size)
scipy.sparse.save_npz(f'{file}.tr_nn_v_cord.rw_beta_csrmat.npz', tr_nn.rw_beta)
tr_nn.rw_beta=scipy.sparse.load_npz(f'{file}.tr_nn_v_cord.rw_beta_csrmat.npz')

olga_beta = OlgaModel(chain_folder = "human_T_beta", recomb_type = "VDJ")
tr_nn.clone_df['pgen_cdr3_b_aa'] = olga_beta.compute_aa_cdr3_pgens(CDR3_seq=tr_nn.clone_df.cdr3_b_aa.to_list())

# Load in VDJdb file and reorganise
fp_vdjdb = 'VDJdb_human_MHCI_MHCII_070524.tsv'
vdjdb = pd.read_csv(fp_vdjdb, sep ='\t')
vdjdb_to_tcrdist = {
    'CDR3':'cdr3_b_aa',
    'V':'v_b_gene',
    'J':'j_b_gene',
    'Score':'score',
    'Species':'species',
    'MHC A':'mhc_a',
    'MHC B':'mhc_b',
    'MHC class':'mhc_class',
    'Epitope':'epitope',
    'Epitope species':'epitope_species'
    }

vdjdb = vdjdb.rename(columns = vdjdb_to_tcrdist)[vdjdb_to_tcrdist.values()]
vdjdb = vdjdb.query('score > 0').reset_index(drop = True)

tr_vdjdb = TCRrep(
    cell_df = vdjdb,
    organism = 'human',
    chains = ['beta'],
    deduplicate = False,
    compute_distances = False)

tr_nn.cpus = 4
chunk_size = get_safe_chunk(tr_nn.clone_df.shape[0],
                            tr_vdjdb.clone_df.shape[0])

tr_nn.compute_sparse_rect_distances(df2 = tr_vdjdb.clone_df,
                                    radius = 72,
                                    chunk_size = chunk_size)

df_join_nn = join_by_dist(
    how = 'left',
    csrmat = tr_nn.rw_beta,
    left_df = tr_nn.clone_df,
    right_df = tr_vdjdb.clone_df,
    left_cols = tr_nn.clone_df.columns.to_list(),
    right_cols = tr_vdjdb.clone_df.columns.to_list(),
    left_suffix = '',
    right_suffix = '_vdjdb',
    max_n = 5,
    radius = 48)

df_join_nn.to_csv('all_hp_bal_vs_vdjdb.tsv', sep = "\t")
 
