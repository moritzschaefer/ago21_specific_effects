import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn2

d = 1000

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
mutant = 'Ago12'

df = pd.read_excel('../../data/TableS3_ATAC-seq_DA_genes.xlsx', skiprows=2)

ago21_specifics = pd.read_excel('../../data/TableS1_RNA-seq.xlsx', skiprows=2, index_col=0, sheet_name='Ago2&1_KO specific DEGs')

gene_l2fcs = df[['Associate GeneID', 'log2FC (DEG)']].drop_duplicates().set_index('Associate GeneID')

df_down = df.query('`DA Status` == "down"')
df_up = df.query('`DA Status` == "up"')
pval = 0.05

for df, label, ax in zip((df_up, df_down), ('up', 'down'), axes):
    promoter_close_genes = df[(df['DA region<->Gene Distance'] < 0) & (df['DA region<->Gene Distance'] > -d) & (df['Pvalue (DA)'] <= pval)]
    da_set = set(promoter_close_genes['Associate GeneID'].dropna().drop_duplicates().values)
    deg_set = set(ago21_specifics.query(f'`Status` == "{label.upper()}"').index)


    venn2((deg_set, da_set), (f'Ago2&1_KO DEGs {label}', f'Genes with {label} chrom. acc.'), ax=ax)

# now plot overlap with ago-specific genes
plt.suptitle('Overlap of Ago2&1 specific DEGs with genes with diff. access. at promoter')
fig.savefig('atac_da_genes_overlap.svg')

