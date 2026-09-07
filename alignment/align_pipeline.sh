#! /bin/bash

##############################################
# modified from ARIMA GENOMICS MAPPING PIPELINE 02/08/2019 #
##############################################

#Below find the commands used to map HiC data.

#Replace the variables at the top with the correct paths for the locations of files/programs on your system.

#This bash script will map one paired end HiC dataset (read1 & read2 fastqs). Feel to modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# inputs #
##########################################
#direcotry containing paired files. Must end _1.fq.gz and _2.fq.gz
IN_DIR=$1
#out directory
OUTPUT=$2
#label for naming files
LABEL=$3
#hg38 or hg19
GENOME=$4
#statement if slurm files are run 1 yes
RUNSLURM=$5
#reuiqred scripts relative to this
SOURCE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

EXCLUDE="--exclude=na05,na46,na48,na50,na53,na42,na49,na06"

#these are the modules loaded in bash scripts
#module load compbio/bwa
#module load compbio/samtools
#module load compbio/picard
#
#for large datasets picard markduplicates seems to require alot of memory. 
PICARD_HIGH_MEM="java -jar -Xms80G -Xmx80G -XX:-UseGCOverheadLimit ./apps/picard-2.27.1/picard.jar"
#reference files

REF=.reference/BWA_${GENOME}_index/$GENOME.fa
FAIDX=.reference/$GENOME.fa.fai

#SRA='basename_of_fastq_files'
#LABEL='overall_exp_name'
#BWA='/software/bwa/bwa-0.7.12/bwa'
#SAMTOOLS='/software/samtools/samtools-1.3.1/samtools'
#IN_DIR='/path/to/gzipped/fastq/files'
#PREFIX='bwa_index_name'

FILTER="$SOURCE_DIR/filter_five_end.pl"
COMBINER="$SOURCE_DIR/two_read_bam_combiner.pl"
STATS="$SOURCE_DIR/get_stats.pl"
#BAM2JUICER="$SOURCE_DIR/bam2juicer.sh"
BAM2JUICER="$SOURCE_DIR/bam2juicer.sh"
JUICER_TOOLS="./juicer/juicer_tools/juicer_tools_1.14.08.jar"

#outputs. Final output is on paired_nodups
RAW_DIR=$OUTPUT/raw_bams/$LABEL
FILT_DIR=$OUTPUT/filtered_bams/$LABEL
TMP_DIR=$OUTPUT/tmp/$LABEL
PAIR_DIR=$OUTPUT/paired/$LABEL
REP_DIR=$OUTPUT/paired_nodups/$LABEL
SCRIPT_DIR=$OUTPUT/scripts/$LABEL
HIC_DIR=$OUTPUT/hic-files
#PICARD='/software/picard/picard-2.6.0/build/libs/picard.jar'
MAPQ_FILTER=30
CPU=12

echo "### Step 0: Check output directories exist & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $SCRIPT_DIR ] || mkdir -p $SCRIPT_DIR
[ -d $HIC_DIR ] || mkdir -p $HIC_DIR
#[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

JOB_ID_STRING=""

for SRA in $IN_DIR/*_1.fq.gz
do
	SRA=${SRA#$IN_DIR/}
	SRA=${SRA%_1.fq.gz}
	
	
	echo "### Step 1.A: FASTQ to BAM (1st) and filter"
	ALIGNER_BATCH1=$SCRIPT_DIR/${SRA}_batch_align1.sh
	echo -n > $ALIGNER_BATCH1
	echo "#!/bin/bash -l
#SBATCH -J ${SRA}_hic_align1
#SBATCH -o ${SRA}_hic_align1_%j.txt
#SBATCH -e ${SRA}_hic_align1_err_%j.txt
#SBATCH -t 7-0
#SBATCH --cpus-per-task=12
#SBATCH -p normal
#SBATCH --mem-per-cpu=12000
#SBATCH --oversubscribe

module load compbio/bwa
module load compbio/samtools


bwa mem -t 12 $REF $IN_DIR/${SRA}_1.fq.gz | samtools view -@ 12 -Sb > $RAW_DIR/${SRA}_1.bam

samtools view -h $RAW_DIR/${SRA}_1.bam | perl $FILTER | samtools view -Sb > $FILT_DIR/${SRA}_1.bam
	" >> $ALIGNER_BATCH1

	# --exclude=na01,na02,na05,na46,na48,na50,na53,na42,na49,na06
if [ $RUNSLURM == 1 ] 
then
	jid1=$(sbatch  $ALIGNER_BATCH1 $EXCLUDE | cut -f 4 -d' ')
	echo $jid1
fi


	echo "### Step 1.B: FASTQ to BAM (1st) and filter"
	ALIGNER_BATCH2=$SCRIPT_DIR/${SRA}_batch_align2.sh
	echo -n > $ALIGNER_BATCH2
	echo "#!/bin/bash -l
#SBATCH -J ${SRA}_hic_align2
#SBATCH -o ${SRA}_hic_align2_%j.txt
#SBATCH -e ${SRA}_hic_align2_err_%j.txt
#SBATCH -t 7-0
#SBATCH --cpus-per-task=12
#SBATCH -p normal
#SBATCH --mem-per-cpu=12000
#SBATCH --oversubscribe

module load compbio/bwa
module load compbio/samtools


bwa mem -t 12 $REF $IN_DIR/${SRA}_2.fq.gz | samtools view -@ 12 -Sb > $RAW_DIR/${SRA}_2.bam

samtools view -h $RAW_DIR/${SRA}_2.bam | perl $FILTER | samtools view -Sb > $FILT_DIR/${SRA}_2.bam
	" >> $ALIGNER_BATCH2
	
if [ $RUNSLURM == 1 ] 
then
	jid2=$(sbatch $ALIGNER_BATCH2 $EXCLUDE | cut -f 4 -d' ')
	echo $jid2
fi


	
	echo "### Step 2: Pair reads & mapping quality filter"
	COMBINER_BATCH=$SCRIPT_DIR/${SRA}_batch_combine.sh
	echo -n > $COMBINER_BATCH
	echo "#!/bin/bash -l
#SBATCH -J ${SRA}_hic_combine
#SBATCH -o ${SRA}_hic_combine_%j.txt
#SBATCH -e ${SRA}_hic_combine_err_%j.txt
#SBATCH -t 7-0
#SBATCH --cpus-per-task=12
#SBATCH -p normal
#SBATCH --mem-per-cpu=8000
#SBATCH --oversubscribe



module load compbio/samtools
module load compbio/picard



perl $COMBINER $FILT_DIR/${SRA}_1.bam $FILT_DIR/${SRA}_2.bam samtools $MAPQ_FILTER | samtools view -bS -t $FAIDX | samtools sort -@ 12 -m 2G -o $TMP_DIR/$SRA.bam
picard AddOrReplaceReadGroups INPUT=$TMP_DIR/$SRA.bam OUTPUT=$PAIR_DIR/${SRA}.bam ID=$SRA LB=$LABEL SM=$LABEL PL=ILLUMINA PU=none
	" >> $COMBINER_BATCH

	if [ $RUNSLURM == 1 ] 
	then
		jid=$(sbatch --dependency=afterok:$jid1,$jid2 $COMBINER_BATCH $EXCLUDE | cut -f 4 -d' ')
		JOB_ID_STRING=$JOB_ID_STRING,$jid
	fi
done



if [ $RUNSLURM == 1 ] 
then
	JOB_ID_STRING=${JOB_ID_STRING#,}
	echo $JOB_ID_STRING
fi
	
	
echo "### Step 3: merge replicates and generate final outputs"
MERGE_BATCH=$SCRIPT_DIR/batch_merge.sh
echo -n > $MERGE_BATCH
echo "#!/bin/bash -l
#SBATCH -J ${LABEL}_hic_merge
#SBATCH -o ${LABEL}_hic_merge_%j.txt
#SBATCH -e ${LABEL}_hic_merge_err_%j.txt
#SBATCH -t 7-0
#SBATCH --cpus-per-task=12
#SBATCH -p normal
#SBATCH --mem-per-cpu=8000
#SBATCH --oversubscribe



module load compbio/samtools
module load compbio/picard

#INPUTS_TECH_REPS=()

#for file in $PAIR_DIR/*.bam; do INPUTS_TECH_REPS+=(\"INPUT=\$file\");done


#picard MergeSamFiles \$(echo \${INPUTS_TECH_REPS[@]}) OUTPUT=$TMP_DIR/$LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

#$PICARD_HIGH_MEM MarkDuplicates INPUT=$TMP_DIR/$LABEL.bam OUTPUT=$REP_DIR/$LABEL.bam METRICS_FILE=$REP_DIR/metrics.$LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

#samtools index $REP_DIR/$LABEL.bam
#perl $STATS $REP_DIR/$LABEL.bam > $REP_DIR/$LABEL.bam.stats

#samtools sort -n $REP_DIR/$LABEL.bam -m 2G -@ 12 -o $TMP_DIR/$LABEL.name.sort.bam

#sh $BAM2JUICER $TMP_DIR/$LABEL.name.sort.bam $REP_DIR/${LABEL}_merged_nodups.txt $TMP_DIR

java -jar -Xms64G -Xmx64G $JUICER_TOOLS pre $REP_DIR/${LABEL}_merged_nodups.txt.gz $HIC_DIR/$LABEL.hic hg38
">>$MERGE_BATCH

if [ $RUNSLURM == 1 ] 
then
	jid=$(sbatch --dependency=afterok:$JOB_ID_STRING $MERGE_BATCH $EXCLUDE)
	echo $jid
fi

