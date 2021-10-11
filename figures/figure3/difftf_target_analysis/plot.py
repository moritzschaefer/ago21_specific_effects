import matplotlib.pyplot as plt
import pandas as pd
import pyensembl
import seaborn as sns

ensembl_release = pyensembl.EnsemblRelease(98, pyensembl.species.mouse)
def geneid_to_name(gene_id):
    '''
    Use like this:
    df['symbol'] = df['geneid'].apply(geneid_to_name)
    '''
    try:
        return ensembl_release.gene_name_of_gene_id(gene_id)
    except (ValueError, KeyError) as e:
        if 'Run pyensembl install' in e.args[0]:
            raise
        else:
            return None
gid2n = geneid_to_name

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

ago21_specifics = pd.read_excel('../../../data/TableS1_RNA-seq.xlsx', skiprows=2, index_col=0, sheet_name='Ago2&1_KO specific DEGs').index
ago21_deg = pd.read_excel('../../../data/TableS1_RNA-seq.xlsx', skiprows=2, index_col=[0, 1], header=[0, 1])['Ago2&1'].droplevel(1)
PVALUE = 0.05
MAX_DISTANCE = 1000

# First, promoter-based target identification

# see main.sh for generation of this file
tfbs_df = pd.read_csv('WTvsAgo12.allMotifs_genes_tss.sorted.merged.bed')
tfbs_df['gene_id'] = tfbs_df['target_id'].apply(lambda v: ensembl_release.transcript_by_id(v).gene_id)
tfbs_df['gene_name'] = tfbs_df['gene_id'].apply(gid2n)
tfbs_df = tfbs_df[(tfbs_df['peak_match_distance'] < MAX_DISTANCE) & (tfbs_df['pval'] < PVALUE)]
# only show the 5 TFs that showed very high statistical significance
TFS = ['ERR2', 'MYC', 'CTCF', 'REST', 'KLF4']
tfbs_df = tfbs_df[tfbs_df.TF.isin(TFS)]

fig, axes = plt.subplots(1, len(TFS), figsize=(len(TFS) * 5.5, 5.5), sharex=True, sharey=True)
full_tss_sets = {}
dfs = {}

tfbs_df['ago21_deg_l2fc'] = ago21_deg.reindex(tfbs_df['gene_id'])['log2FoldChange'].fillna(0).values
tfbs_df['ago21_deg_padj'] = ago21_deg.reindex(tfbs_df['gene_id'])['padj'].fillna(1.0).values
for ax, (tf, tf_df) in zip(axes, tfbs_df.groupby('TF')):
    per_gene = tf_df.groupby('gene_id')[['l2FC', 'ago21_deg_l2fc', 'ago21_deg_padj']].mean()
    per_gene['Ago21-specific'] = pd.Categorical(['no'] * len(per_gene), categories=['no', 'yes'], ordered=True)
    per_gene.loc[per_gene.index.intersection(ago21_specifics), 'Ago21-specific'] = 'yes'
    explained = ((per_gene['l2FC'] > 0) == (per_gene['ago21_deg_l2fc'] > 0)) & (per_gene['ago21_deg_padj'] < 0.05)
    per_gene['Gene name'] = per_gene.index.map(gid2n)
    per_gene['explained'] = explained
    per_gene['source'] = 'promoter'
    dfs[tf] = [per_gene]
    pos_explainable_genes = (explained & (per_gene['l2FC'] > 0))
    neg_explainable_genes = (explained & (per_gene['l2FC'] < 0))
    sns.scatterplot(data=per_gene.sort_values('Ago21-specific'), y='l2FC', x='ago21_deg_l2fc', hue='Ago21-specific', ax=ax, palette={'yes': '#Ff8c00', 'no': '#Cecece'})
    ax.text(ax.get_xlim()[1] * 0.6, ax.get_ylim()[1] * 0.6, (pos_explainable_genes & (per_gene['Ago21-specific'] == 'no')).sum())
    ax.text(ax.get_xlim()[1] * 0.6, ax.get_ylim()[1] * 0.5, (pos_explainable_genes & (per_gene['Ago21-specific'] == 'yes')).sum(), color='#Ff8c00')

    ax.text(ax.get_xlim()[0] * 0.7, ax.get_ylim()[0] * 0.7, (neg_explainable_genes & (per_gene['Ago21-specific'] == 'no')).sum())
    ax.text(ax.get_xlim()[0] * 0.7, ax.get_ylim()[0] * 0.8, (neg_explainable_genes & (per_gene['Ago21-specific'] == 'yes')).sum(), color='#Ff8c00')
    ax.set_ylabel('TFBS DA (l2fc)')
    label = tf
    ax.set_title(label)
    full_tss_sets[label] = set(pos_explainable_genes.index[pos_explainable_genes]) | set(neg_explainable_genes.index[neg_explainable_genes])
_ = fig.suptitle('Gene-centric TFBS analysis. Number of "explained" genes are indicated in the corners (DA and DEG are stat. sign.).')
fig.savefig('difftf_targets_deg_integration.svg')

# Second, enhancer-based target identification

MAX_DISTANCE = 0
fig, axes = plt.subplots(1, len(TFS), figsize=(len(TFS) * 5.5, 5), sharex=True, sharey=True)
table_df = pd.read_csv('WTvsAgo12.allMotifs_enhancers_merged.merged.bed')

if table_df['target_id'].iloc[0].startswith('ENSMUST'):
    table_df['gene_id'] = table_df['target_id'].apply(lambda v: ensembl_release.transcript_by_id(v).gene_id)
    table_df['gene_name'] = table_df['gene_id'].apply(gid2n)
else:
    table_df['gene_name'] = table_df['target_id']
    table_df['gene_id'] = table_df['gene_name'].apply(gn2id)

# TODO do we have negative distances here? that would be bad, because we only want to allow overlaps with enhancer regions (i.e. i need to add .abs())
table_df = table_df[(table_df['peak_match_distance'] <= MAX_DISTANCE) & (table_df['pval'] < PVALUE)]
TFS = ['ERR2', 'MYC', 'CTCF', 'REST', 'KLF4']
table_df = table_df[table_df.TF.isin(TFS)]
table_df['ago21_deg_l2fc'] = ago21_deg.reindex(table_df['gene_id'])['log2FoldChange'].fillna(0).values
table_df['ago21_deg_padj'] = ago21_deg.reindex(table_df['gene_id'])['padj'].fillna(1.0).values

for i, (ax, (tf, tf_df)) in enumerate(zip(axes, table_df.groupby('TF'))):
    label = tf
    ax.set_title(label)

    per_gene = tf_df.groupby('gene_id')[['l2FC', 'ago21_deg_l2fc', 'ago21_deg_padj']].mean()
    per_gene['Ago21-specific'] = pd.Categorical(['no'] * len(per_gene), categories=['no', 'yes'], ordered=True)
    per_gene.loc[per_gene.index.intersection(ago21_specifics), 'Ago21-specific'] = 'yes'
    per_gene['Gene name'] = per_gene.index.map(gid2n)
    explained = ((per_gene['l2FC'] > 0) == (per_gene['ago21_deg_l2fc'] > 0)) & (per_gene['ago21_deg_padj'] < 0.05)
    per_gene['explained'] = explained
    per_gene['source'] = 'enhancer'
    dfs[tf].append(per_gene)
    pos_explainable_genes = (explained & (per_gene['l2FC'] > 0))
    neg_explainable_genes = (explained & (per_gene['l2FC'] < 0))
    xx = sns.scatterplot(data=per_gene.sort_values('Ago21-specific'), y='l2FC', x='ago21_deg_l2fc', hue='Ago21-specific', ax=ax, palette={'yes': '#Ff8c00', 'no': '#Cecece'})
    if i == 0:
        xx.legend(loc='lower right')
    else:
        xx.legend().remove()

    ax.text(ax.get_xlim()[1] * 0.7, ax.get_ylim()[1] * 0.7, (pos_explainable_genes & (per_gene['Ago21-specific'] == 'no')).sum())
    ax.text(ax.get_xlim()[1] * 0.7, ax.get_ylim()[1] * 0.5, (pos_explainable_genes & (per_gene['Ago21-specific'] == 'yes')).sum(), color='#Ff8c00')

    ax.text(ax.get_xlim()[0] * 0.7, ax.get_ylim()[0] * 0.65, (neg_explainable_genes & (per_gene['Ago21-specific'] == 'no')).sum())
    ax.text(ax.get_xlim()[0] * 0.7, ax.get_ylim()[0] * 0.95, (neg_explainable_genes & (per_gene['Ago21-specific'] == 'yes')).sum(), color='#Ff8c00')
    if i == 0:
        ax.set_ylabel('TFBS DA (l2fc)')
    sns.despine(ax=ax)

_ = fig.suptitle('Gene-centric TFBS analysis. Number of "explained" genes are indicated in the corners (DA and DEG are stat. sign.).')
fig.savefig('difftf_targets_deg_integration_published_annotations.svg')
