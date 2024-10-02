#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 18:58:48 2024

@author: wezi
"""

import pandas as pd
from tcrdist.repertoire import TCRrep
from tcrdist.public import TCRpublic

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

##
## MAKE PUBLICITY REPORT AND OUTPUT AS HTML
##

# <tp> TCRpublic class for reporting publicities 
tp = TCRpublic(
    tcrrep = tr_bulk, 
    output_html_name = "all_hp_bal_publicity.html")
# set to True, if we want a universal radius
tp.fixed_radius = True
# must then specify maximum distance for finding similar TCRs
tp.radius = 18
# set criteria for being quasi-public
tp.query_str = 'nsubject > 1'
# Add additional columns to be summarized in the report
tp.kargs_member_summ['addl_cols'] = ['subject']

# by calling, .report() an html report is made
public = tp.report()

##
## ABOVE SECTION IS TO MAKE PUBLICITY REPORT
##