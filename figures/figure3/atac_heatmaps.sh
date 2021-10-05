#!/usr/bin/env bash

# transcript regions
computeMatrix scale-regions -p 4 -R ../../data/genes.bed -a 3000 -b 3000 --regionBodyLength 8000 --skipZeros -S ../../pipelines/atac_seq/subsampled/{WT_1,WT_2,Ago12_1,Ago12_2}.filtered.bam.bw -o gene_regions_matrix.mat.gz
plotHeatmap -m gene_regions_mat.mat.gz -o gene_regions.png

# promoters
computeMatrix reference-point --referencePoint TSS -p 4 -R ../../data/genes.bed -a 3000 -b 3000 --skipZeros -S ../../{WT_1,WT_2,Ago12_1,Ago12_2}.filtered.bam.bw -o promoters.mat.gz
plotProfile -m promoters.mat.gz --perGroup -o profile_promoters.png &

# enhancers were compiled from Table S3 and S4 from Gonz\'alez-Ram\'irez, Mar, Ballar\'e, Cecilia, Mugianesi, F., Beringer, M., Santanach, A., Blanco, E., & Di Croce, L., Differential contribution to gene expression prediction of histone modifications at enhancers or promoters, PLOS Computational Biology, 17(9), 1009368 (2021).  http://dx.doi.org/10.1371/journal.pcbi.1009368

computeMatrix reference-point --referencePoint TSS -p 4 -R ../../data/genes.bed -a 3000 -b 3000 --skipZeros -S ../../{WT_1,WT_2,Ago12_1,Ago12_2}.filtered.bam.bw -o enhancers.mat.gz
plotProfile -m enhancers.mat.gz --perGroup -o profile_enhancers.png &

wait
