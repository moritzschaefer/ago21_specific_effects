#!/bin/bash
# requires snakePipes
createIndices --local -o $HOME/data/snakepipes/GRCm38_98 --genome ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz GRCm38_98 -j 25
