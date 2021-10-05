import pandas as pd

df = pd.read_csv('encode_h3k27me3_ensembl.wig', sep='\t')

scaled = (df / 3)
scaled=scaled.dropna()
scaled = scaled.reset_index()
scaled['level_1'] = scaled['level_1'].astype(int)
scaled['level_2'] = scaled['level_2'].astype(int)
scaled['level_0'] = scaled['level_0'].apply(lambda v: f'chr{v}')
scaled.to_csv('encode_h3k27me3_ensembl.scaled.wig', header=False, sep='\t', index=False)
