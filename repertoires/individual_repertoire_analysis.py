#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 13:49:43 2024

@author: wezi
"""
import sys
import pandas as pd
import numpy as np
import scipy.sparse
from tcrdist.breadth import get_safe_chunk
from scipy.stats import poisson
from tcrdist.repertoire import TCRrep
from tcrdist.join import join_by_dist
from tcrdist.public import _neighbors_sparse_fixed_radius
from tcrdist.background import sample_britanova
from tcrdist.background import get_stratified_gene_usage_frequency
from statsmodels.stats.multitest import multipletests
from tcrdist.pgen import OlgaModel

CPUS = 4

# load up the repertoire specified in command line argument
# (they're the .tsv files in immunarch folder)
file = sys.argv[1]

df_bulk = pd.read_csv(file, sep='\t')

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

tr_bulk.weights_b = {'cdr3_b_aa': 6, 'pmhc_b_aa': 1, 'cdr2_b_aa': 1, 'cdr1_b_aa': 1}
chunk_size = get_safe_chunk(tr_bulk.clone_df.shape[0], tr_bulk.clone_df.shape[0])
tr_bulk.compute_sparse_rect_distances(radius = 72, chunk_size = chunk_size)
scipy.sparse.save_npz(f'{file}.tr_bulk.rw_beta_csrmat.npz', tr_bulk.rw_beta)

tr_bulk.rw_beta=scipy.sparse.load_npz(f'{file}.tr_bulk.rw_beta_csrmat.npz')

# Looking up within-group neigbours. TINKER AROUND WITH RADIUS HERE
tr_bulk.clone_df['nn'] = _neighbors_sparse_fixed_radius(
    csrmat = tr_bulk.rw_beta,
    radius=48)

# Count neighbours
tr_bulk.clone_df['k_nn'] = [len(x)-1 for x in tr_bulk.clone_df['nn']]

# Uncommenting this would define dataframe of clones with five or more
# within-group neighbours
# df_nn = tr_bulk.clone_df.iloc[tr_bulk.clone_df.query('k_nn > 4').index,].reset_index(drop = True)

# load up cord blood sample
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

tr_nn.clone_df['nn_cord'] = _neighbors_sparse_fixed_radius(
    csrmat = tr_nn.rw_beta,
    radius = 48)

tr_nn.clone_df['k_nn_cord'] = [len(x) for x in tr_nn.clone_df['nn_cord']]

bulk_n = tr_bulk.clone_df.shape[0]
Q = 1
tr_nn.clone_df['lambda'] = (tr_nn.clone_df['k_nn_cord']+1)/tr_cord.clone_df.shape[0]
tr_nn.clone_df['poisson'] = tr_nn.clone_df.apply(lambda x : 1- poisson.cdf(x['k_nn'], Q*x['lambda']*bulk_n, loc = 0), axis = 1)
tr_nn.clone_df['poisson_fdr'] = multipletests(tr_nn.clone_df['poisson'], method = "fdr_bh")[1]
tr_nn.clone_df['poisson_holm'] = multipletests(tr_nn.clone_df['poisson'], method = "holm")[1]
tr_nn.clone_df['poisson_bonferroni'] = multipletests(tr_nn.clone_df['poisson'], method = "bonferroni")[1]

np.sum(tr_nn.clone_df['poisson_fdr'] < 0.001)

tr_nn.clone_df[['v_b_gene', 'cdr3_b_aa', 'k_nn', 'k_nn_cord', 'poisson_fdr']].sort_values('poisson_fdr')

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

output_filename = file + '_vs_vdjdb.tsv'
df_join_nn.to_csv(output_filename, sep = "\t")
 
