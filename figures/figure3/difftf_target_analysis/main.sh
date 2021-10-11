#!/usr/bin/env bash

# enhancer data is from Gonz\'alez-Ram\'irez, Mar, Ballar\'e, Cecilia, Mugianesi, F., Beringer, M., Santanach, A., Blanco, E., & Di Croce, L., Differential contribution to gene expression prediction of histone modifications at enhancers or promoters, PLOS Computational Biology, 17(9), 1009368 (2021).  http://dx.doi.org/10.1371/journal.pcbi.1009368
BEDFILES="../../../data/genes_tss.sorted.bed journal.pcbi.1009368.s020_TableS8_active_enhancers_refined.bed journal.pcbi.1009368.s021_TableS9_poised_enhancers_refined.bed"

# sort bedfiles
for bedfile in $BEDFILES ../difftf/diffTF_repo/output/FINAL_OUTPUT/extension100/WTvsAgo12.allMotifs.bed; do
    sort -k1,1 -k2,2n $bedfile > $bedfile.sorted &
done
wait

for bedfile in $BEDFILES; do
    bedtools closest -d -a ../difftf/diffTF_repo/output/FINAL_OUTPUT/extension100/WTvsAgo12.allMotifs.bed.sorted \
             -b $bedfile.sorted > WTvsAgo12.allMotifs_$(basename $bedfile) &
done
wait

# combine active and poised enhancer data
cp WTvsAgo12.allMotifs_journal.pcbi.1009368.s020_TableS8_active_enhancers_refined.bed WTvsAgo12.allMotifs_enhancers_merged.bed
tail +2 WTvsAgo12.allMotifs_journal.pcbi.1009368.s021_TableS9_poised_enhancers_refined.bed >> WTvsAgo12.allMotifs_enhancers_merged.bed

python associate_peaks_genes.py
python plot.py
