'''
Plot peak diffexp. vs peak size

To run this pipeline you need to run the atac_seq pipeline first!
'''
import itertools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import binned_statistic

ATAC_MUTANTS = ['Ago12', 'Ago1', 'Ago2']
FDR_FILTER=0.05
LOG2FC_FILTER=0.3

COVERAGES = [f'../../pipelines/atac_seq/coverages/{mutant}_coverage_short_Genrich.bed' for mutant in ATAC_MUTANTS]
DIFF_BINDINGS = [f'../../pipelines/atac_seq/CSAW_Genrich_{mutant}_diff_atac/DiffBinding_scores.txt' for mutant in ATAC_MUTANTS]
SCALE_FACTORS = [f'../../pipelines/atac_seq/manual_scale_factors_{mutant.lower()}.tsv' for mutant in ATAC_MUTANTS]

palette = {'up': 'salmon', 'down': 'lightskyblue', 'insignificant': 'darkgray'}

fig, axes = plt.subplots(2, 3, figsize=(15, 11), sharex=True, sharey=True)
for ax, mutant, cov_file, diff_binding_file, scale_factor_file in zip(
        itertools.chain(*axes),
        ATAC_MUTANTS,
        COVERAGES,
        DIFF_BINDINGS,
        SCALE_FACTORS):

    data_cols = [mutant + '_1', mutant + '_2', 'WT_1', 'WT_2']
    region_counts = pd.read_csv(
        cov_file,
        header=None,
        names=['chr', 'start', 'end', 'name',
               'bindingScore', 'ignore1'] + data_cols,
        sep='\t').set_index('name')

    # either normalize by both, or just normalize by the factors. first by both:
    scale_factors = pd.read_csv(scale_factor_file, sep='\t', index_col=0)
    region_counts.iloc[:, -4:] *= scale_factors.values.reshape(1, -1)

    # binding scores
    df = pd.read_csv(diff_binding_file, sep='\t').set_index('name')
    mean_log2_counts = (region_counts.iloc[:, -4:] + 1) \
        .transform(np.log2).sum(axis=1)/4

    # this showed very high correlation with 'best.logFC'
    manual_log2fc = np.log2(region_counts.iloc[:, -4:-2].sum(axis=1) \
                            / region_counts.iloc[:, -2:].sum(axis=1))
    df['peak log2FC'] = manual_log2fc
    df['mean log2 peak counts'] = mean_log2_counts
    df['direction'] = manual_log2fc.apply(lambda v: 'up' if v > 0 else 'down')
    df.loc[df.FDR >= FDR_FILTER, 'direction'] = 'insignificant'
    df.loc[df['peak log2FC'].abs() < LOG2FC_FILTER, 'direction'] = 'insignificant'

    # deleting these is justified since they mostly correspond to tiny peaks
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    # sort by direction so that insignificant is drawn first
    df['direction'] = pd.Categorical(df['direction'],
                                     categories=['insignificant',
                                                 'up', 'down'],
                                     ordered=True)
    df = df.sort_values('direction')
    df = df.join(region_counts[['chr', 'start', 'end']])

    sns.scatterplot(data=df, y='peak log2FC', x='mean log2 peak counts',
                    hue='direction', palette=palette, ax=ax, s=11,
                    rasterized=True,
                    linewidth=0, alpha=1)

    mean_log2fc_data = binned_statistic(df['mean log2 peak counts'],
                                        values=df['peak log2FC'].to_numpy(),
                                        bins=100,
                                        range=(1,
                                               np.quantile(mean_log2_counts,
                                                           0.9997)))

    ax.plot(mean_log2fc_data.bin_edges[:-1],
            np.nan_to_num(mean_log2fc_data.statistic),
            color='dimgray', linewidth=3)
    ax.set_title(mutant)
    sns.despine(ax=ax)
    ax.grid()
    ax.legend().remove()

    ax.text(7,  3.5, f'n={(df["direction"] == "up").sum()}', color=palette['up'])
    ax.text(7, -3.5, f'n={(df["direction"] == "down").sum()}', color=palette['down'],
            verticalalignment='top')
    ax.text(13, 1, f'n={len(df)}', color=palette['insignificant'])


plt.tight_layout()
fig.savefig('atac_ma_plots.svg')
