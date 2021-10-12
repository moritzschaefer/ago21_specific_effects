import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn2

d = 1000

fig, ax = plt.subplots(1, 1, figsize=(7, 6))
mutant = 'Ago12'

df = pd.read_excel('../../data/TableS3_ATAC-seq_DA_genes.xlsx', skiprows=2)
gene_l2fcs = df[['Associate GeneID', 'log2FC (DEG)']].drop_duplicates().set_index('Associate GeneID')

df_down = df.query('`DA Status` == "down"')
df_up = df.query('`DA Status` == "up"')
ax.axvline(0, color='black')
ax.axhline(0.5, color='black')

pval = 0.05
for df, label, color, spec_color in zip((df_down, df_up), ('dec.', 'inc.'), ('#87cefa', '#fa8072'), ('#98f5ff', '#Ff7f00')):
    promoter_close_genes = df[(df['DA region<->Gene Distance'] < 0) & (df['DA region<->Gene Distance'] > -d) & (df['Pvalue (DA)'] <= pval)]
    l2fc = gene_l2fcs.loc[promoter_close_genes['Associate GeneID'].dropna().drop_duplicates().values]
    l2fc = l2fc.where(l2fc != 0.0)  # in our table, we described missing log2FCs as 0.0. Here we don't want to show them.
    print(color, spec_color)
    sns.ecdfplot(l2fc, label=f'{len(l2fc)} genes with {label} chrom. access. at promoter (<{d}bp, DA Pval<{pval})', ax=ax, color=color)

    # plot specific-filtered
    l2fc = gene_l2fcs.loc[promoter_close_genes.query('`Ago2&1_KO Specific` == "yes"')['Associate GeneID'].dropna().drop_duplicates().values]
    sns.ecdfplot(l2fc, label=f'{len(l2fc)} {label}-DA genes, ago21-specific-filtered', ax=ax, color=spec_color)

ax.set_xlim([-1, 1])
ax.legend()
ax.grid(False)
sns.despine()
ax.set_title(mutant)

plt.suptitle('Diff. expression of genes with differential chrom. access. at promoters')
fig.savefig('atac_da_regions_deg_cdf.svg')

# now plot overlap with ago-specific genes
