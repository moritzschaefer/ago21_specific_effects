#!/bin/bash

# requires snakePipes

mkdir coverages 2> /dev/null

DNA-mapping  -i /mnt/cclab_nas/rawdata/ATAC-seq/fastq -o . --fastqc --mapq 3 -j 32 --dedup GRCm38_98 --local --trim --trimmer trimgalore --trimmerOptions "--nextera --paired"

MUTANTS="Ago12 Ago1 Ago2"

# first align peaks to genes
for mutant in $MUTANTS; do
    cp manual_scale_factors_${mutant}.tsv manual_scale_factors.tsv
    for peak_caller in Genrich; do # HMMRATAC MACS2
        ATAC-seq -d . -c config.yaml --local --sampleSheet ${mutant}_diff_atac.tsv --peakCaller $peak_caller --snakemakeOptions '\--verbose --printshellcmds' -j 25 GRCm38_98
    done
    # if coverage for peak-regions are required
    bedtools multicov -bed CSAW_MACS2_${mutant}_diff_atac/DiffBinding_allregions.bed -bams short_bams/${mutant}_{1,2}.short.cleaned.bam short_bams/WT_{1,2}.short.cleaned.bam > coverages/${mutant}_coverage_short_${peak_caller}.bed
done

# subsampling factors were calculate as described in README.txt
mkdir subsampled
cp filtered_bam/Ago12_1.filtered.bam subsampled/
samtools view -s 0.617 -o subsampled/WT_2.filtered.bam filtered_bam/WT_2.filtered.bam & 
samtools view -s 0.633 -o subsampled/WT_1.filtered.bam filtered_bam/WT_1.filtered.bam &
samtools view -s 0.83333 -o subsampled/Ago12_2.filtered.bam filtered_bam/Ago12_2.filtered.bam & 
samtools view -s 0.706  -o subsampled/Ago2_2.filtered.bam filtered_bam/Ago2_2.filtered.bam & 
samtools view -s 0.77  -o subsampled/Ago2_1.filtered.bam filtered_bam/Ago2_1.filtered.bam & 
samtools view -s 0.696 -o subsampled/Ago1_2.filtered.bam filtered_bam/Ago1_2.filtered.bam &
cd subsampled
for f in *.bam; do
    bamCoverage -b $f -o $f.bw
done
