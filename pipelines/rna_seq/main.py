#!/usr/bin/env python
# requires snakePipes

import glob
import os
import subprocess
from pathlib import Path

import pandas as pd

try:
    os.mkdir('./data')
    # generate symlinks for all rna-seq files into one folder
    # files.csv was copied from moritzschaefers rna-seq-star-deseq2 pipeline
    sample_df = pd.read_csv('files.csv')
    for _, row in sample_df.iterrows():
        for i in range(1, 3):
            os.symlink(row[f'fq{i}'], f'data/{row["sample"]}_R{i}.fastq.gz')
except FileExistsError:
    pass

for sheet in glob.glob('sampleSheets/*.tsv'):
    subprocess.check_call([
    'mRNA-seq', '-i', './data', '-o', '.', '--local', '--trim', '--trimmer', 'trimgalore', '--trimmerOptions', '--illumina --paired', '-j', '32', '--mode', 'alignment,alignment-free,deepTools_qc', '--sampleSheet', sheet, '--fastqc', 'GRCm38_98'])
