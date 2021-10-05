#!/bin/bash
set -ex


ALLSAMPLES=("20190607.A-1_WT_Max_K27_R1" "20190607.A-2_2i_WT_Max_K27_R1" "20190607.A-3_Ago1_KO1_K27_R1" "20190607.A-4_Ago2_KO1_K27_R1" "20190607.A-5_Input_WT_Max_R1" "20190607.A-6_Input_2i_WT_Max_R1" "20190607.A-7_Input_Ago1_KO1_R1" "20190607.A-8_Input_Ago2_KO1_R1" "20190607.A-9_Dgcr8_KO1_K27_R1" "20190607.A-10_Drosha_KO1_K27_R1" "20190607.A-11_Dicer_KO1_K27_R1" "20190607.A-12_Ago21_KO1_K27_R1" "20190607.A-13_Input_Dgcr8_KO1_R1" "20190607.A-14_Input_Drosha_KO1_R1" "20190607.A-15_Input_Dicer_KO1_R1" "20190607.A-16_Input_Ago21_KO1_R1" "20190607.A-17_WT_Dani_K27_R1" "20190607.A-18_2i_WT_Dani_K27_R1" "20190607.A-19_Ago1_KO2_K27_R1" "20190607.A-20_Ago2_KO2_K27_R1" "20190607.A-21_Input_WT_Dani_R1" "20190607.A-22_Input_2i_WT_Dani_R1" "20190607.A-23_Input_Ago1_KO2_R1" "20190607.A-24_Input_Ago2_KO2_R1" "20190607.A-25_Dgcr8_KO2_K27_R1" "20190607.A-26_Drosha_KO2_K27_R1" "20190607.A-27_Dicer_KO2_K27_R1" "20190607.A-28_Ago21_KO2_K27_R1" "20190607.A-29_Input_Dgcr8_KO2_R1" "20190607.A-30_Input_Drosha_KO2_R1" "20190607.A-31_Input_Dicer_KO2_R1" "20190607.A-32_Input_Ago21_KO2_R1")

IP_SAMPLES=("20190607.A-1_WT_Max_K27_R1" "20190607.A-2_2i_WT_Max_K27_R1" "20190607.A-3_Ago1_KO1_K27_R1" "20190607.A-4_Ago2_KO1_K27_R1" "20190607.A-9_Dgcr8_KO1_K27_R1" "20190607.A-10_Drosha_KO1_K27_R1" "20190607.A-11_Dicer_KO1_K27_R1" "20190607.A-12_Ago21_KO1_K27_R1" "20190607.A-17_WT_Dani_K27_R1" "20190607.A-18_2i_WT_Dani_K27_R1" "20190607.A-19_Ago1_KO2_K27_R1" "20190607.A-20_Ago2_KO2_K27_R1" "20190607.A-25_Dgcr8_KO2_K27_R1" "20190607.A-26_Drosha_KO2_K27_R1" "20190607.A-27_Dicer_KO2_K27_R1" "20190607.A-28_Ago21_KO2_K27_R1")

IP_SAMPLES_2=("20190607.A-1_WT_Max_K27_R1" "20190607.A-2_2i_WT_Max_K27_R1" "20190607.A-3_Ago1_KO1_K27_R1" "20190607.A-4_Ago2_KO1_K27_R1" "20190607.A-12_Ago21_KO1_K27_R1")

# set directories
BASEDIR="/path/to/basedir"
FASTQ="$BASEDIR/fastq"
TRIMMED_FASTQ="$BASEDIR/fastq/trimmed"
BAMFILES="$BASEDIR/bam/"
BIGWIGFILES="$BASEDIR/big/"
BEDGRAPHFILES="$BASEDIR/bedgraph/"
HEATMAPSFILES="$BASEDIR/heatmaps/"
MATRIXFILES="$BASEDIR/matrix/"
SUMMARYFILES="$BASEDIR/summary_files"
GENE_ID="$BASEDIR/gene_id"
TRIMMOATIC="/path/to/trimmomatic.jar"
ADAPTER_FILE="$BASEDIR/TruSeq3-SE.fa"

#set the Genome Reference index
INDEX="$BASEDIR/../mm10/mm10"
DM_INDEX="$BASEDIR/../dm6/genome"

#Trimming
trim() {
  pushd $FASTQ
  for file in $(ls grep .fastq.gz)
  do
      java -jar $TRIMMOMATIC SE -threads 10 -phred33 $file \
        $TRIMMED_FASTQ/${file%.fastq.gz}_trimmed.fastq.gz \
        ILLUMINACLIP:$ADAPTER_FILE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
  done
  popd
}

#Alignment
align() {
  echo "Aligning samples..."
  for SAMPLE in "${ALLSAMPLES[@]}"; do
    echo $SAMPLE
    bowtie2 -p 8 -x $1 -U $TRIMMED_FASTQ/${SAMPLE}_trimmed.fastq.gz 2> $BASEDIR/logs/${SAMPLE}_${2}_trimmed_align.log | samtools view -q 10 -h -b -F0x104 - > $BAMFILES/${SAMPLE}_${2}_trimmed.bam # https://www.biostars.org/p/271006/ F0x100 = secondary alignment = when aligned two times take one away, F0x4 = segment unmapped = take away unmapped segment
   
    # Adapted from markdup documentation (http://www.htslib.org/doc/samtools.html) instead of using rmdup = because this command is obsolete/
    # fixmate, namesorted and mate scoring was left out because we have no paired sequencing
    samtools sort -@ 8 -T tmpSort -m 6G -o $BAMFILES/${SAMPLE}_${2}_sorted.bam $BAMFILES/${SAMPLE}_${2}_trimmed.bam
   
  
    # Finally mark duplicates
    samtools markdup -r -s $BAMFILES/${SAMPLE}_${2}_sorted.bam $BAMFILES/${SAMPLE}_${2}_trimmed_final.bam
    samtools index $BAMFILES/${SAMPLE}_${2}_trimmed_final.bam
    
    rm $BAMFILES/${SAMPLE}_${2}_sorted.bam 
    
  done
}
#### Index/DM_index corresponds to numer 1 and mm/dm corresponds to number 2 --> that way it is possible to only write one alignment command and not one for mm and one for dm

## compare ip over input with bamCompare ans spike in
compare(){
    echo "bam compare"
    IP_INDICES=(1 2 3 4 9 10 11 12 17 18 19 20 25 26 27 28)
  
  for IP_INDEX in "${IP_INDICES[@]}"; do  
    IP_NAME=$(echo "${BAMFILES}/20190607.A-${IP_INDEX}_"*"_mm_trimmed_final.bam")
    INPUT_NAME=$(echo "$BAMFILES/20190607.A-$(($IP_INDEX + 4))_Input_*_mm_trimmed_final.bam")
    SPIKEIN_NAME=$(echo "${BAMFILES}/20190607.A-${IP_INDEX}_"*"_dm_trimmed_final.bam")
    
    OUTPUT_NAME=${IP_NAME##*/}
    OUTPUT_NAME=${OUTPUT_NAME%.*} 
    
    SAMPLE_DM_COUNT=$(samtools view -c $SPIKEIN_NAME) 
    INPUT_COUNT=$(samtools view -c $INPUT_NAME)
    INPUT_ADJUSTED=$(bc <<< "scale=3; $INPUT_COUNT/507.358") # factor calculated as described in readme

    bamCompare -p 8 -b1 $IP_NAME -b2 $INPUT_NAME -o "$BIGWIGFILES/${OUTPUT_NAME}_trimmed_normalized.bigwig" -of bigwig --scaleFactors "1:$(bc <<< "scale=3; $INPUT_ADJUSTED/$SAMPLE_DM_COUNT")"
  done
}
####--blackListFileName  --> no blacklist used in this script

compute_matrix() {
  for SAMPLE in "${IP_SAMPLES_2[@]}"; do ## GTF file for start of gene regions downloaded form gencode https://www.gencodegenes.org/mouse/
    INPUT_NAME=$(echo "${BIGWIGFILES}/${SAMPLE}_mm_trimmed_final_trimmed_normalized.bigwig") #corresponds to input used for script not CHIP INPUT!
    OUTPUT_NAME=${INPUT_NAME##*/}
    OUTPUT_NAME=${OUTPUT_NAME%.*}
  
    computeMatrix reference-point -S $INPUT_NAME -R $BASEDIR/gencode.vM23.annotation.gtf -a 5000 -b 5000 -o "$MATRIXFILES/${OUTPUT_NAME}_test.matrix" -p 8 #--outFileNameMatrix "$MATRIXFILES/${OUTPUT_NAME}_test.tab" #--transcriptID "exon" # in output name next time include folder e.g. "$MATRIXFILES/${OUTPUT_NAME}.matrix"
  done
}

#computeMatrixOperations cbind -m "$MATRIXFILES/20190607.A-1_WT_Max_K27_R1_mm_trimmed_final_trimmed_normalized.matrix" "$MATRIXFILES/20190607.A-12_Ago21_KO1_K27_R1_mm_trimmed_final_trimmed_normalized.matrix" "$MATRIXFILES/20190607.A-3_Ago1_KO1_K27_R1_mm_trimmed_final_trimmed_normalized.matrix" "$MATRIXFILES/20190607.A-4_Ago2_KO1_K27_R1_mm_trimmed_final_trimmed_normalized.matrix"  -o "$MATRIXFILES/IP_samples_2_combined_test.matrix" ##to combine matrices

plot_Heatmap(){
  echo "plot Heatmaps"
 for SAMPLE in "${IP_SAMPLES_2[@]}"; do 
 
   INPUT_NAME=$(echo "${MATRIXFILES}/${SAMPLE}_mm_trimmed_final_trimmed_normalized.matrix") #corresponds to input used for script not CHIP INPUT!
   OUTPUT_NAME=${INPUT_NAME##*/}
   OUTPUT_NAME=${OUTPUT_NAME%.*}
 
   plotHeatmap -m $INPUT_NAME -o "$HEATMAPSFILES/${OUTPUT_NAME}_test.png" --colorMap Purples --kmeans 4  #--zMin 0 --zMax 2.5 # #--clusterUsingSamples 1 #--outFileSortedRegions $GENE_ID/${OUTPUT_NAME}_test.bed  #  --kmeans 4 
 done
}

## Quality check of Chiped samples
bam_summary_files() {
  echo "bam_summary_files"
   IP_INDICES=(1 2 3 4 9 10 11 12 17 18 19 20 25 26 27 28)
  
  for IP_INDEX in "${IP_INDICES[@]}"; do  
    IP_NAME=$(echo "${BAMFILES}/20190607.A-${IP_INDEX}_"*"_mm_final.bam")
    #INPUT_NAME=$(echo "$BAMFILES/20190607.A-$(($IP_INDEX + 4))_Input_*_mm_final.bam")
    
    OUTPUT_NAME=${IP_NAME##*/}
    OUTPUT_NAME=${OUTPUT_NAME%.*}
  
    multiBamSummary bins --bamfiles $BAMFILES/20190607.A-1_WT_Max_K27_R1_mm_final.bam $BAMFILES/20190607.A-2_2i_WT_Max_K27_R1_mm_final.bam $BAMFILES/20190607.A-3_Ago1_KO1_K27_R1_mm_final.bam $BAMFILES/20190607.A-4_Ago2_KO1_K27_R1_mm_final.bam $BAMFILES/20190607.A-9_Dgcr8_KO1_K27_R1_mm_final.bam $BAMFILES/20190607.A-10_Drosha_KO1_K27_R1_mm_final.bam $BAMFILES/20190607.A-11_Dicer_KO1_K27_R1_mm_final.bam $BAMFILES/20190607.A-12_Ago21_KO1_K27_R1_mm_final.bam $BAMFILES/20190607.A-17_WT_Dani_K27_R1_mm_final.bam $BAMFILES/20190607.A-18_2i_WT_Dani_K27_R1_mm_final.bam $BAMFILES/20190607.A-19_Ago1_KO2_K27_R1_mm_final.bam $BAMFILES/20190607.A-20_Ago2_KO2_K27_R1_mm_final.bam $BAMFILES/20190607.A-25_Dgcr8_KO2_K27_R1_mm_final.bam $BAMFILES/20190607.A-26_Drosha_KO2_K27_R1_mm_final.bam $BAMFILES/20190607.A-27_Dicer_KO2_K27_R1_mm_final.bam $BAMFILES/20190607.A-28_Ago21_KO2_K27_R1_mm_final.bam -o $SUMMARYFILES/test_IP.npz -p 6
  
  done
}

## multiBamSummary bins --bamfiles $IP_NAME -o $SUMMARYFILES/test_IP.npz
## plotPCA -in $SUMMARYFILES/test_IP.npz -o $SUMMARYFILES/PCA_IP.png -T "PCA of IP"


# Ip indices of Input = Ip Index + 4, thats why we can use the numbers instead of the full name
peak_call () { # one argument: if 0 then no spikein, if 1 spikein
  echo "peak calling "
  IP_INDICES=(1 2 3 4 9 10 11 12 17 18 19 20 25 26 27 28)
  
  for IP_INDEX in "${IP_INDICES[@]}"; do  
    INPUT_NAME=$(echo "${BAMFILES}/20190607.A-${IP_INDEX}_"*"_mm_trimmed_final.bam")
    OUTPUT_NAME=$(echo "${BAMFILES}/20190607.A-${IP_INDEX}_"*"_mm_trimmed_final.bam")
    
    OUTPUT_NAME=${INPUT_NAME##*/}
    OUTPUT_NAME=${OUTPUT_NAME%.*} 
 
    macs2 callpeak -t $INPUT_NAME -c $BAMFILES/20190607.A-$(($IP_INDEX + 4))_Input_*_mm_trimmed_final.bam --broad -g mm --broad-cutoff 0.1  --outdir $BASEDIR/peaks -n ${OUTPUT_NAME} -f BAM 
  done  
}

#merging bigwig files: (need to have installed: ucsc-bigwigmerge, ucsc-bedgraphtobigwig, ucsc-fetchchromsizes)

#bigWigMerge "$BIGWIGFILES/20190607.A-1_WT_Max_K27_R1_mm_trimmed_final_trimmed_normalized.bigwig" "$BIGWIGFILES/20190607.A-17_WT_Dani_K27_R1_mm_trimmed_final_trimmed_normalized.bigwig" "$BEDGRAPHFILES/test_WT_Max_Dani_merged.bedGraph"
#fetchChromSizes mm10 > mm10.chrom.sizes
#sort -k1,1 -k2,2n "BEDGRAPHFILES/test_WT_Max_Dani_merged.bedGraph" > "BEDGRAPHFILES/sorted_merged_WT_MAX.bedGraph
#bedGraphToBigWig "$BEDGRAPHFILES/test_WT_Max_Dani_merged.bedGraph" mm10.chrom.sizes "BIGWIGFILES/test.bw"

#other possibility use bigwigcompare from deeptools:
#bigwigCompare -b1 "$BIGWIGFILES/20190607.A-1_WT_Max_K27_R1_mm_trimmed_final_trimmed_normalized.bigwig" -b2 "$BIGWIGFILES/20190607.A-17_WT_Dani_K27_R1_mm_trimmed_final_trimmed_normalized.bigwig" -o test.bigwig -of bigwig -p 8

main() {
  trim
  align $INDEX "mm"
  align $DM_INDEX "dm"
  minimum_drosophila_count
  coverage
  coverage_spike_in
  compare
  compute_matrix
  plot_Heatmap
  plot_Fingerprint
  bam_summary_files
  peak_call
}

main











