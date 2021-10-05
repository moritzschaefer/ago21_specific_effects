#!/usr/bin/env bash

# TODO Need to actually compute and use genes_tss.bed instead of genes.bed
# TODO maybe needs to be sorted
# enhancer data is from Gonz\'alez-Ram\'irez, Mar, Ballar\'e, Cecilia, Mugianesi, F., Beringer, M., Santanach, A., Blanco, E., & Di Croce, L., Differential contribution to gene expression prediction of histone modifications at enhancers or promoters, PLOS Computational Biology, 17(9), 1009368 (2021).  http://dx.doi.org/10.1371/journal.pcbi.1009368
for bedfile in ../../data/genes.bed journal.pcbi.1009368.s020_TableS8_active_enhancers_refined.bed journal.pcbi.1009368.s021_TableS9_poised_enhancers_refined.bed; do
    bedtools closest -d -a ../difftf/diffTF_repo/output_extension100/WTvsAgo12.allMotifs.bed \
             -b $bedfile > WTvsAgo12.allMotifs_$(basename $bedfile)
done

cp WTvsAgo12.allMotifs_journal.pcbi.1009368.s020_TableS8_active_enhancers_refined.bed WTvsAgo12.allMotifs_enhancers_merged.bed
tail +2 WTvsAgo12.allMotifs_journal.pcbi.1009368.s021_TableS9_poised_enhancers_refined.bed >> WTvsAgo12.allMotifs_enhancers_merged.bed

python associate_peaks_genes.py
python plot.py
