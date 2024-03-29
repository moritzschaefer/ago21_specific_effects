#!/bin/bash

# wget https://www.encodeproject.org/files/ENCFF296RAS/@@download/ENCFF296RAS.bigWig
bigWigToWig ENCFF296RAS.bigWig encode_h3k27me3_ensembl.wig

# scale to our WT sample
python scale_encode_wig.py

# convert back to bigwig (chrom.sizes is provided here for GRCm38.98 but can also be downloaded from ENSEMBL or UCSC and generated manually)
wigToBigWig encode_h3k27me3_ensembl.scaled.wig files/chrom.sizes.sorted encode_h3k27me3_ensembl.scaled.bw

# bigwig files (other than the encode one) are generated by the h3k27me3-pipeline and can also be obtained from GEO
computeMatrix reference-point --referencePoint TSS -S encode_h3k27me3_ensembl.scaled.bw 20190607.A-1_WT_Max_K27_R1_mm_trimmed_final_trimmed_normalized.bigwig 20190607.A-12_Ago21_KO1_K27_R1_mm_trimmed_final_trimmed_normalized.bigwig -R ../../data/genes.bed --beforeRegionStartLength 5000  --afterRegionStartLength 5000 -p 4 --skipZeros -o tss_with_encode.normalized.mat.gz

# Generate clusters using our WT sample. clustered transcripts are stored as bed-file
plotHeatmap --matrixFile tss_with_encode.normalized.mat.gz -out tss_with_encode.normalized.png --outFileSortedRegions tss_with_encode.normalized.bed --clusterUsingSamples 2 --kmeans 4
