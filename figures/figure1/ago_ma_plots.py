import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plotlib import maplot, volcanoplot

_parentdir = pathlib.Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(_parentdir))


mrna_df = pd.read_excel('../../data/TableS1_RNA-seq.xlsx', index_col=[0, 1], header=[0, 1], skiprows=2).droplevel(level=1, axis=0)
ago21_specifics = pd.read_excel('../../data/TableS1_RNA-seq.xlsx', skiprows=2, index_col=0, sheet_name='Ago2&1_KO specific DEGs').index
mutants = ('Ago1', 'Ago2', 'Ago2&1')

fig, axes = plt.subplots(2, 3, figsize=(12, 12))

for ax, mutant in zip(np.transpose(axes), mutants):
    maplot(mrna_df[mutant], mean='tpm_expression', ax=ax[0])
    plotdf = mrna_df[mutant].copy()
    plotdf['-log10(adjusted pvalue)'] = -np.log10(mrna_df[mutant]['padj'])
    volcanoplot(plotdf, ax=ax[1], show_n=True)
    if mutant == 'Ago2&1':
        volcanoplot(plotdf.loc[ago21_specifics], ax=ax[1], show_n=True, palette={'up': 'orange', 'down': 'green', 'unchanged': 'black'})

    ax[0].set_title(mutant)

fig.savefig('ago_ma_volcano_plots.svg')
