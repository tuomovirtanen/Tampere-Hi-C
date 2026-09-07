#!/bin/bash
#SBATCH --job-name=count_reads_mt
#SBATCH --output=logs/count_reads_mt.%A_%a.out
#SBATCH --error=logs/count_reads_mt.%A_%a.err
#SBATCH --partition=large
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-44    # 45 samples total (0-based index)

# Load modules
module load SAMtools
module load BEDTools

# Set number of threads
THREADS=4

files=(
BPH_337.bam
BPH_456.bam
BPH_651.bam
BPH_652.bam
BPH_656.bam
BPH_659.bam
BPH_671.bam
BPH_677.bam
BPH_688.bam
BPH_689.bam
BPH_701.bam
CRPC_261.bam
CRPC_278.bam
CRPC_305.bam
CRPC_348.bam
CRPC_435.bam
CRPC_489.bam
CRPC_539.bam
CRPC_541.bam
CRPC_542.bam
CRPC_543.bam
CRPC_697.bam
PC_12517.bam
PC_14670.bam
PC_15194.bam
PC_15420.bam
PC_15760.bam
PC_17163.bam
PC_17447.bam
PC_18307.bam
PC_19403.bam
PC_20873.bam
PC_22392.bam
PC_22603.bam
PC_24173.bam
PC_470.bam
PC_4786.bam
PC_4980.bam
PC_6174.bam
PC_6342.bam
PC_6488.bam
PC_7875.bam
PC_8131.bam
PC_8438.bam
PC_9324.bam
)

BAM=${files[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename "$BAM" .bam)

BAM_PATH="./Tampere_prostate_HiC/atac/bams"
OUTPATH="./Hi-C/ATAC-peak-analysis/results/counts"
BEDPATH="./Hi-C/ATAC-peak-analysis/results/bins_100000_masked.sorted.bed"
PEAKPATH="./Hi-C/ATAC-peak-analysis/results/peaks.sorted.bed"
#CHROMSIZE_FILE="./Hi-C/ATAC-new/scripts/hg38.chrom.sizes.sorted.txt"
echo "Processing $BAM"

#the bam files are soretd this makes sure the chromsizes are in the same order, The non canonical contigs seemed to be in weird ordrer, bu the canonical ones were correct?
CHROMSIZE_FILE="$OUTPATH/${SAMPLE}_hg38.from_bam.chrom.sizes"
samtools view -H ./Tampere_prostate_HiC/atac/bams/BPH_337.bam | grep '^@SQ' | awk -F'\t' '{OFS="\t"; gsub("SN:", "", $2); gsub("LN:", "", $3); print $2, $3}' > $CHROMSIZE_FILE

# Run the command
#1280 = 256 + 1024 i.e. exclude secondary alignments and duplicated reads
#modify these e.g. if you want to keep duplicates
#MAPQ 30
#-f 2 properly paired reads only
samtools view -b -q 30 -f 2 -F 1280 -@ $THREADS "$BAM_PATH/$BAM" | bedtools coverage -sorted -g $CHROMSIZE_FILE -a $BEDPATH -b - -counts | gzip > "$OUTPATH/${SAMPLE}_bin_counts.bed.gz"
samtools view -b -q 30 -f 2 -F 1280 -@ $THREADS "$BAM_PATH/$BAM" | bedtools coverage -sorted -g $CHROMSIZE_FILE -a $PEAKPATH -b - -counts | gzip > "$OUTPATH/${SAMPLE}_peak_counts.bed.gz"

# 
echo "Done with $SAMPLE"

