import matplotlib.pyplot as plt
import pandas as pd
import pyensembl
import requests
from matplotlib_venn import venn3
from pybedtools import BedTool

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


def http_get(url, fn):
    with open(fn, 'wb') as f:
        f.write(requests.get(url, allow_redirects=True).content)

# load DEG sets
ago21_deg = pd.read_excel('../../data/TableS1_RNA-seq_V3.xlsx', skiprows=2, index_col=[0, 1]).droplevel(1)['Ago2&1']
ago21_up_degs = set(ago21_deg.query('log2FoldChange > 0.5 & padj < 0.05').index)
ago21_down_degs = set(ago21_deg.query('log2FoldChange < -0.5 & padj < 0.05').index)
specific_genes = pd.read_excel('../../data/TableS1_RNA-seq_V3.xlsx', skiprows=2, index_col=0, sheet_name='Ago2&1_KO specific DEGs').index
specific_up = ago21_up_degs & set(specific_genes)
specific_down = ago21_down_degs & set(specific_genes)

# load processed motif bindings (for TSS matches)
# TODO get and gunzip
motifs = pd.read_csv('../figure3/difftf/diffTF_repo/output/FINAL_OUTPUT/extension100/WTvsAgo12.allMotifs.tsv.gz', gunzip=True, sep='\t')
motifs = motifs[['chr', 'MSS', 'MES', 'permutation', 'TF', 'TFBSID', 'strand', 'peakID',
       'l2FC', 'limma_avgExpr', 'limma_B', 'limma_t_stat', 'DESeq_ldcSE',
       'DESeq_stat', 'DESeq_baseMean', 'pval', 'pval_adj']]

# TODO generate (README.md)
tss_bed = BedTool('../../data/genes_tss.sorted.bed')

# load enhancer-gene bindings
enhancer_sites = pd.concat([
    pd.read_csv(f, sep='\t', header=False, names=['chr_e', 'start_e', 'end_e', 'chr_p', 'start_p', 'end_p', 'gene']) for f in
    ['../figure3/difftf_target_analysis/journal.pcbi.1009368.s020_TableS8_active_enhancers_refined.bed'
'../figure3/difftf_target_analysis/journal.pcbi.1009368.s021_TableS9_poised_enhancers_refined.bed']
])
enhancer_sites['chr_e'] = enhancer_sites['chr_e'].apply(lambda v: v.replace('chr', ''))
enhancer_bed = BedTool.from_dataframe(enhancer_sites)

# plot CTCF
# data from: Nora, E. P., Goloborodko, A., Valton, A., Gibcus, J. H., Uebersohn, A., Abdennur, N., Dekker, J., …, Targeted degradation of ctcf decouples local insulation of chromosome domains from genomic compartmentalization, Cell, 169(5), 930–944–22 (2017).  http://dx.doi.org/10.1016/j.cell.2017.05.004

http_get('https://ars.els-cdn.com/content/image/1-s2.0-S0092867417305317-mmc4.xlsx', 'ctcf_chip.xlsx')

ctcf_df = pd.read_excel('ctcf_chip.xlsx', skiprows=6, header=None)[[0, 1, 2]].dropna()
ctcf_df.columns = ['chr', 'start', 'end']
ctcf_df['start'] = ctcf_df['start'].astype(int)
ctcf_df['end'] = ctcf_df['end'].astype(int)
ctcf_df['chr'] = ctcf_df['chr'].astype(str).str.replace('chr', '')

ctcf_bed = BedTool.from_dataframe(ctcf_df)

ctcf_motifs = motifs[motifs['TF'] == 'CTCF']
da_ctcf_motifs = ctcf_motifs[ctcf_motifs['pval'] < 0.05]

ctcf_tss_matches = BedTool.from_dataframe(ctcf_df).window(tss_bed, w=1000).to_dataframe(header=None)
ctcf_tss_matches['gene_id'] = ctcf_tss_matches[6].apply(lambda tid: ensembl_release.transcript_by_id(tid).gene_id)

# additionally, only look into peaks that were not detected by diffTF (because we previously studied diffTF targets already)
other_inter = ctcf_bed.intersect(BedTool.from_dataframe(ctcf_motifs), u=True).to_dataframe()
other_inter['chrom'] = other_inter['chrom'].astype(str)

ctcf_enhancer_matches = BedTool.from_dataframe(ctcf_df).window(enhancer_bed, w=0).to_dataframe()
ctcf_enhancer_matches['gene_id'] = ctcf_enhancer_matches.iloc[:, -1].apply(gn2id)

fig, axes = plt.subplots(1, 2, figsize=(15, 5))
venn3((specific_up, set(ctcf_tss_matches.gene_id) & ago21_up_degs, set(ctcf_enhancer_matches.gene_id) & ago21_up_degs), ('Ago21 Specific up genes', 'DEGs with TSS-proximal peaks', 'DEGs with enhancer-proximal peaks'), ax=axes[0])
venn3((specific_down, set(ctcf_tss_matches.gene_id) & ago21_down_degs, set(ctcf_enhancer_matches.gene_id) & ago21_down_degs), ('Ago21 Specific down genes', 'DEGs with TSS-proximal peaks', 'DEGs with enhancer-proximal peaks'), ax=axes[1])
axes[0].set_title('Overlap with upregulated DEGs')
axes[1].set_title('Overlap with downregulated DEGs')
fig.suptitle('Overlap of Ctcf-ChIP-seq-peak associated genes with Ago21-specific DEGs')
fig.savefig('ctcf_chip_overlap_ago21degs.svg')

# plot KLF4
http_get('https://static-content.springer.com/esm/art%3A10.1038%2Fs41556-019-0390-6/MediaObjects/41556_2019_390_MOESM3_ESM.xlsx', 'klf4_chip.xlsx')
# filter transient because there, the peaks are very weak in mESC
klf4_df = pd.read_excel('klf4_chip.xlsx', skiprows=1).query('cluster != "Transient"')
klf4_df['chr'] = klf4_df['chr'].astype(str).str.replace('chr', '')

klf4_bed = BedTool.from_dataframe(klf4_df)
klf4_motifs = motifs[motifs['TF'] == 'KLF4']
da_klf4_motifs = klf4_motifs[klf4_motifs['pval'] < 0.05]

klf4_tss_matches = BedTool.from_dataframe(klf4_df).window(tss_bed, w=1000).to_dataframe(header=None)
klf4_tss_matches['gene_id'] = klf4_tss_matches[7].apply(lambda tid: ensembl_release.transcript_by_id(tid).gene_id)
klf4_enhancer_matches = BedTool.from_dataframe(klf4_df).window(enhancer_bed, w=0).to_dataframe()
klf4_enhancer_matches['gene_id'] = klf4_enhancer_matches.iloc[:, -1].apply(gn2id)

fig, axes = plt.subplots(1, 2, figsize=(15, 5))
venn3((specific_up, set(klf4_tss_matches.gene_id) & ago21_up_degs, set(klf4_enhancer_matches.gene_id) & ago21_up_degs), ('Ago21 Specific up genes', 'DEGs with TSS-proximal peaks', 'DEGs with enhancer-proximal peaks'), ax=axes[0])
venn3((specific_down, set(klf4_tss_matches.gene_id) & ago21_down_degs, set(klf4_enhancer_matches.gene_id) & ago21_down_degs), ('Ago21 Specific down genes', 'DEGs with TSS-proximal peaks', 'DEGs with enhancer-proximal peaks'), ax=axes[1])
axes[0].set_title('Overlap with upregulated DEGs')
axes[1].set_title('Overlap with downregulated DEGs')
fig.suptitle('Overlap of Klf4-ChIP-seq-peak associated genes with Ago21-specific DEGs')
fig.savefig('klf4_chip_overlap_ago21degs.svg')
