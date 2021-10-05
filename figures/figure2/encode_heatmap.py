
import glob
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyensembl
import seaborn as sns

ensembl_release = pyensembl.EnsemblRelease(98, pyensembl.species.mouse)

PROMOTER_DISTANCE = 5000

print('To run this file, please download the histone mark ChIP-seq files from ENCODE as indicated in the Materials&Methods section of our publication')

### Prepare gene/transcript locations
# genes.bed contains all transcripts as provided by ENSEMBL GRCm38 v98
df = pd.read_csv('../../data/genes.bed', sep='\t', header=None)
df[0] = df[0].astype(str)
df[12] = df[3].apply(ensembl_release.gene_name_of_transcript_id)
CHROMOSOMES = [str(v) for v in range(1, 21)] + ['X', 'Y']

### add promoter-region (5000 bp) to each transcript region
drop_indices = []
for i, row in df.iterrows():
    if row[5] == '+':
        df.loc[i, 1] = max(0, df.loc[i, 1] - PROMOTER_DISTANCE)
    else:
        df.loc[i, 2] += PROMOTER_DISTANCE

    if row[0] not in CHROMOSOMES:
        drop_indices.append(i)
df.drop(drop_indices, inplace=True)
df.set_index(12, inplace=True)

### convert gene groups into bed files
for f in glob.glob('files/*.csv'):
    if 'result' in f:
        continue

    group = pd.read_csv(f)
    df.loc[df.index.intersection(group['gene_name'])].to_csv(f.replace('.csv', '.bed'), index=False, sep='\t', header=False)

### Get histone mark levels at gene/transcript regions
for f in glob.glob('files/*.bed'):
    subprocess.check_call(['multiBigwigSummary', 'BED-file', '--outFileName', f.replace('.bed', '.npz'), '--BED', f, '--bwfiles'] + list(glob.glob('files/*.bw')))

### Take mean over each group for each histone mark
means = {}

for group in glob.glob('files/*.npz'):
    name = Path(group).stem

    data = np.load(group)
    df = pd.DataFrame(data.get('matrix'), columns=data.get('labels'))
    means[name] = df.fillna(0).mean()
df = pd.DataFrame(means).T

df.columns = df.columns.map(lambda v: v.replace('.bw', ''))  # .replace('_zscore', '')
df.index = df.index.map(lambda v: v.replace('genes', 'All genes').replace('expressed_All genes', 'expressed genes').replace('_DEGs', '').replace('RNA_seq_', '').replace('_volcano', '').replace('_', ' ').title())

df = df.loc[[
    'Ago21 Specific Up',
    'Ago21 Specific Down',
    'Mirna Targets Schaefer2021',
    'Expressed Genes',
    'All Genes'
]]
df = df[['H3K9me3', 'H3K27me3', 'H3K4me1', 'H3K27ac', 'H3K9ac', 'H3K36me3', 'H3K4me3']]

# normalize zscore
for col in df.columns:
    df[col] = (df[col] - df[col].mean())/df[col].std(ddof=0)

fig, ax = plt.subplots(figsize=(6, len(df)/1.4))
sns.heatmap(df, center=0, cmap='vlag', cbar_kws={'label': 'column-wise z-score'})
fig.savefig('encode_histone_heatmap.png')
