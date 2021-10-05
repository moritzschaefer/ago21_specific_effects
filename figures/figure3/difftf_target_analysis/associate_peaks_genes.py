# coding: utf-8
import glob

import pandas as pd
import pyensembl

df = pd.read_csv('WTvsAgo12.allMotifs.tsv', sep='\t')
df.set_index('TFBSID', inplace=True)

# this code serves for different kinds of bed files. Here we only use one for now
for f in ['WTvsAgo12.allMotifs_genes.bed', 'WTvsAgo12.allMotifs_enhancers_merged.bed']:
    matches = pd.read_csv(f, sep='\t', header=None)
    if matches.iloc[0, 9].startswith('ENSMUST'):
        matches.rename(columns={9: 'target_id', 18: 'peak_match_distance'}, inplace=True)
    elif matches.iloc[0, 9].startswith('chr'):
        matches.rename(columns={12: 'target_id', 13: 'peak_match_distance'}, inplace=True)
    else:
        matches.rename(columns={9: 'target_id', 10: 'peak_match_distance'}, inplace=True)

    merged = df.merge(matches.set_index(3)[['target_id', 'peak_match_distance']], left_index=True, right_index=True)
    merged['pval_adj'] = merged['pval_adj'].replace('  NA', 1).astype(float)
    merged = merged[merged['peak_match_distance'] < 1000]  # normally 0 should be it, but we can also filter later
    merged.to_csv(f.replace('.bed', '.merged.bed'))
