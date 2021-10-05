import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def maplot(data, mean='baseMean', diff='log2FoldChange', p='padj', p_threshold=0.05, l2fc_threshold=0.5, ax=None, palette={'up': 'red', 'down': 'blue', 'unchanged': 'black'}, label_ss=False):
    data = data.copy()
    data['hue'] = 'unchanged'
    data.loc[(data[p] < p_threshold) & (data[diff] > l2fc_threshold), 'hue'] = 'up'
    data.loc[(data[p] < p_threshold) & (data[diff] < -l2fc_threshold), 'hue'] = 'down'
    data['log2(mean expression + 1)'] = np.log2(data[mean] + 1)
    volcanoplot(data=data, ax=ax, show_n='y', x='log2(mean expression + 1)', y=diff, s=7, palette=palette, label_ss=label_ss)


def volcanoplot(data, x='log2FoldChange', y='-log10(adjusted pvalue)', show_n=False, ax=None, alpha=1.0, hue=None, palette={'up': 'red', 'down': 'blue', 'unchanged': 'black'}, y_threshold=-np.log10(0.05), x_threshold=0.5, s=40, label_ss=False):
    '''
    Given a df, plot a volcano plot.

    Either provide a field for the coloring on your own (hue parameter), or provide a threshold via y_threshold and x_threshold to compute the color

    Default fields: 'log2FoldChange' and '-log10(adjusted pvalue)'
    cutoffs for p-value is 0.05 and for log2FoldChange 0.5
    '''
    data = data.copy()
    if not hue:
        hue = 'hue'
    if hue not in data.columns:
        data[hue] = 'unchanged'
        data.loc[(data[x] > x_threshold) & (data[y] > y_threshold), hue] = 'up'
        data.loc[(data[x] < -x_threshold) & (data[y] > y_threshold), hue] = 'down'

    # sort by hue so that insignificant is drawn first
    data['hue'] = pd.Categorical(data['hue'], categories=['unchanged', 'up', 'down'], ordered=True)
    data = data.sort_values('hue')

    ax = sns.scatterplot(data=data, x=x, y=y, hue=hue, palette=palette, s=s, edgecolor="none",
                         legend=False, alpha=alpha, ax=ax, rasterized=True, linewidth=0)
    if show_n:
        if show_n == 'y':
            x_up, x_down = data[x].max(), data[x].max()
            y_up = min(data[y].max() * 0.8, ax.get_ylim()[1])
            y_down = max(data[y].min() * 0.8, ax.get_ylim()[0])
            ha_up = 'right'
            ha_down = 'right'
        else:
            y_up, y_down = data[y].max(), data[y].max()
            x_up = min(data[x].max() * 0.9, ax.get_xlim()[1])
            x_down = max(data[x].min() * 0.9 , ax.get_xlim()[0])
            ha_up = 'right'
            ha_down = 'left'

        ax.text(x_down, y_down, f'n={(data[hue] == "down").sum()}', color=palette['down'], horizontalalignment=ha_down)
        ax.text(x_up, y_up, f'n={(data[hue] == "up").sum()}', color=palette['up'], horizontalalignment=ha_up)

    if label_ss:
        for index, row in data[~(data[hue] == 'unchanged')].iterrows():
            ax.text(row[x], row[y], index)

    sns.despine()

    return ax
