mkdir ../results
mkdir ../results/counts
mkdir ../results/coverage
mkdir ../results/corrected_peaks

#paths in the scripts are defined relative to this folder.

#When generating the bins. the black lsited peaks are included in the peak file to mask those regions from the bins
#We exclude these peaks in the differential peak analysis later


Rscript generate_bins_and_GC_fractions.r peak_file=../results/peaks.bed outpath=../results

#sort bin bed file to same order as the bam files have been sorted to speed up the counting step
sort -k1,1 -k2,2n ../results/bins_100000_masked.bed > ../results/bins_100000_masked.sorted.bed

#also sort the peaks, the order is not chr1,chr2,chr3 but chr1,chr10,chr11,...
#keep the original indexes for for easy mapping back later
#awk '{print $0, NR}' ../results/peaks.bed | sort -k1,1 -k2,2n > ../results/peaks.sorted.bed
awk 'BEGIN{OFS="\t"}{print $0, NR}' ../results/peaks.bed | sort -k1,1 -k2,2n > ../results/peaks.sorted.bed

#this will calculate the coverage in the bins for all samples using SLURM job array
#modify the paths and samples etc. in this file!!!
sbatch calculate_bin_coverage.sh

#generate peak count matrix from the individual files
Rscript peak_count_matrix.r ../results/counts

#this will generate the log2 ratios in the 100kb windows (or whatever bin size you used)
#check the paths, and the reference pattern in the script before running (here I use BPH samples as reference)
Rscript coverage_log2ratios.r


#this script you need to modify according to your needs
#this will do the differential peak analysis using the peak counts and the log2 ratios form the previous step
#check the paths and the reference pattern in the script before running (here I use BPH samples as reference when calculaing the offsets)
#also you need set the design matrix according to your sample groups in the script
Rscript differential_peak_analysis.r