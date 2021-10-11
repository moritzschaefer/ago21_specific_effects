import matplotlib.pyplot as plt
import pandas as pd
import pyensembl
import requests
from matplotlib_venn import venn2

ensembl_release = pyensembl.EnsemblRelease(98, pyensembl.species.mouse)
def genename_to_id(gene_name):
    '''
    Use like this:
    df['symbol'] = df['geneid'].apply(genename_to_id)
    '''
    try:
        return ensembl_release.gene_ids_of_gene_name(gene_name)[0]
    except (ValueError, KeyError) as e:
        if 'Run pyensembl install' in e.args[0]:
            raise
        else:
            return None


gn2id = genename_to_id

ago21_specifics = pd.read_excel('../../data/TableS1_RNA-seq.xlsx', skiprows=2, sheet_name='Ago2&1_KO specific DEGs')
ups = set(ago21_specifics.query('Status == "UP"')['Gene ID'])
downs = set(ago21_specifics.query('Status == "DOWN"')['Gene ID'])

cluster1 = set(pd.read_excel('../../data/TableS2_H3K27me3_clusters.xlsx', skiprows=2).query('Cluster == "cluster_1"')['Gene ID'])

# From (Asenjo et al. 2020, Table S4)
bivalent_genes = pd.read_excel('aay4768_table_s4.xlsx', index_col=0, sheet_name='HC-Bivalent ChIP-seqs').index  # the other sheet name: 'HC-Bivalent-4sUseq'
bivalent_genes = set(bivalent_genes.map(gn2id))

fig, ax = plt.subplots(1, 1, figsize=(6, 5))
venn2((cluster1, bivalent_genes), ('Cluster_1', 'Bivalent_genes'), ax=ax)
fig.savefig('cluster1_bivalent_overlap.png')

fig, ax = plt.subplots(1, 1, figsize=(6, 5))
venn2((cluster1, ups), ('Cluster_1', 'Ago2&1_KO specific up DEGs'), ax=ax)
fig.savefig('cluster1_ups_overlap.png')

fig, ax = plt.subplots(1, 1, figsize=(6, 5))
venn2((cluster1, downs), ('Cluster_1', 'Ago2&1_KO specific down DEGs'), ax=ax)
fig.savefig('cluster1_downs_overlap.png')
