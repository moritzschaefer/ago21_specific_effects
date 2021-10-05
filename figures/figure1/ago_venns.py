import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn3

groups = {}
mrna_df = pd.read_excel('../../data/TableS1_RNA-seq.xlsx', index_col=[0, 1], header=[0, 1], skiprows=2).droplevel(level=1, axis=0)
mutants = ('Ago1', 'Ago2', 'Ago2&1')
# for mutant in ('Ago1', 'Ago2', 'Ago2&1'):
#     df = pd.read_excel('../../data/TableS1_RNA-seq.xlsx', sheet_name='{mutant}_KO DEGs')
#     groups[mutant] = df

#venn3([set(v) for v in groups.values()], [set(k) for k in groups.keys()])

sets = [set(mrna_df[mutant].query('(log2FoldChange > 0.5 | log2FoldChange < -0.5) & padj < 0.05').index) for mutant in mutants]
fig, ax = plt.subplots(1, 1, figsize=(6, 5))

venn3(sets, mutants)

plt.savefig('ago21_ago1_ago2_overlap.png')
