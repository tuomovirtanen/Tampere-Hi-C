# Material and methods for Tampere PC ATAC project

## Scripts used
The following section lists all the scripts used to process the ATAC-seq FASTQ files up to the quantification of the consensus peaks.

### Trimming FASTQs
FASTQ reads were trimmed using `trim_galore` version 0.6.10.

```bash
trim_galore \
    -j 8 \
    --paired \
    --length 20 \
    -q 20 \
    --output_dir "$OUTDIR" \
    "${PAIRED_FASTQS[@]}"
```

```bash
#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --output=logs/trim_galore/%x_%a.out
#SBATCH --error=logs/trim_galore/%x_%a.err
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-38
#SBATCH --partition=small

set -euo pipefail

module load Trim_Galore

BASE_IN="../data/fastqs"
BASE_OUT="../data/fastqs_trimmed"

# Deterministic sample list
mapfile -t SAMPLES < <(find "$BASE_IN" -mindepth 1 -maxdepth 1 -type d | sort)

SAMPLE_DIR="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"
SAMPLE_NAME="$(basename "$SAMPLE_DIR")"

OUTDIR="${BASE_OUT}/${SAMPLE_NAME}"
mkdir -p "$OUTDIR"

cd "$SAMPLE_DIR"

# Build pairwise-ordered FASTQ list
PAIRED_FASTQS=()

for R1 in *_R1_*.fastq.gz; do
    R2="${R1/_R1_/_R2_}"

    if [[ ! -f "$R2" ]]; then
        echo "ERROR: Missing R2 for $R1"
        exit 1
    fi

    PAIRED_FASTQS+=("$R1" "$R2")
done

if [[ ${#PAIRED_FASTQS[@]} -eq 0 ]]; then
    echo "ERROR: No FASTQs found in $SAMPLE_NAME"
    exit 1
fi

trim_galore \
    -j 8 \
    --paired \
    --length 20 \
    -q 20 \
    --output_dir "$OUTDIR" \
    "${PAIRED_FASTQS[@]}"

```

### Genome alignment
Trimmed reads were aligned to the `hg38` reference genes (analysis set: `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna`) using `bowtie2` version 2.5.1, using `--very-sensitive -X 2000` parameters.

```bash
#!/bin/bash
#SBATCH --job-name=bowtie2_align
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/bowtie2/%x_%a.out
#SBATCH --error=logs/bowtie2/%x_%a.err
#SBATCH --array=1-38
#SBATCH --partition=small

set -euo pipefail

module load Bowtie2
module load SAMtools

BASE_IN="../data/fastqs_trimmed"
BASE_OUT="../data/bams"
INDEX="../references/GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bowtie_index"

mkdir -p "$BASE_OUT"

# Deterministic list of samples
mapfile -t SAMPLES < <(find "$BASE_IN" -mindepth 1 -maxdepth 1 -type d | sort)

SAMPLE_DIR="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"
SAMPLE="$(basename "$SAMPLE_DIR")"

cd "$SAMPLE_DIR"

# Collect trimmed FASTQs and build comma-separated lists
R1_LIST=$(ls *_R1_*_val_1.fq.gz | sort | paste -sd,)
R2_LIST=$(ls *_R2_*_val_2.fq.gz | sort | paste -sd,)

if [[ -z "$R1_LIST" || -z "$R2_LIST" ]]; then
    echo "ERROR: Missing trimmed FASTQs for $SAMPLE"
    exit 1
fi

OUT_BAM="${BASE_OUT}/${SAMPLE}_realigned_sorted.bam"

bowtie2 \
  -x "$INDEX" \
  -1 "$R1_LIST" \
  -2 "$R2_LIST" \
  -p "${SLURM_CPUS_PER_TASK}" \
  --sensitive-local \
  -X 2000 | \
samtools view -@ "${SLURM_CPUS_PER_TASK}" -b | \
samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "$OUT_BAM"

samtools index "$OUT_BAM"
```

### Marking duplicate reads
Duplicate reads were marked (but not removed) using `Picard MarkDuplicates` version VN:3.0.0 using `VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=FALSE`.

```bash
#!/bin/bash
#SBATCH --job-name=markdup
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/picard_markdup/%x_%a.out
#SBATCH --error=logs/picard_markdup/%x_%a.err
#SBATCH --array=1-38
#SBATCH --partition=small

set -euo pipefail

module load picard

BASE_IN="../data/bams"
BASE_OUT="../data/bams_markdup"
TMP_BASE="${BASE_OUT}/tmp"

mkdir -p "$BASE_OUT"
mkdir -p "$TMP_BASE"

# sample list based on BAMs
mapfile -t BAMS < <(find "$BASE_IN" -name "*_realigned_sorted.bam" | sort)

IN_BAM="${BAMS[$((SLURM_ARRAY_TASK_ID - 1))]}"
SAMPLE="$(basename "$IN_BAM" _realigned_sorted.bam)"

echo "SAMPLE=${SAMPLE}"
echo "SAMPLE=${SAMPLE}" >&2

OUT_BAM="${BASE_OUT}/${SAMPLE}_realigned_sorted_markdup.bam"
METRICS="${BASE_OUT}/${SAMPLE}_markdup_metrics.txt"
TMP_DIR="${TMP_BASE}/${SAMPLE}"

mkdir -p "$TMP_DIR"

java -Xmx28g -Djava.io.tmpdir="$TMP_DIR" -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
    INPUT="$IN_BAM" \
    OUTPUT="$OUT_BAM" \
    METRICS_FILE="$METRICS" \
    VALIDATION_STRINGENCY=LENIENT \
    REMOVE_DUPLICATES=FALSE \
    CREATE_INDEX=TRUE \
    TMP_DIR="$TMP_DIR"

```

### Converting the BAM files to BED format
BAM files were converted to BED format to be used by `macs3`. At the same time, removing **duplicate reads**, **secondary alignments**, **supplementary alignments**, and **unmapped reads**, as well as reads wth `MAPQ < 20`.

```bash
#!/bin/bash
#SBATCH --job-name=bam2bed
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3000
#SBATCH --output=logs/bam2bed/%x_%a.out
#SBATCH --error=logs/bam2bed/%x_%a.err
#SBATCH --array=1-38
#SBATCH --partition=small

set -euo pipefail

module load SAMtools
samtools --version

module load BEDTools
bedtools --version

module load all/HTSlib
bgzip --version

BASE_IN="../data/bams_markdup"
BASE_OUT="../data/beds"
TMP_BASE="${BASE_OUT}/tmp"

mkdir -p "$BASE_OUT"
mkdir -p "$TMP_BASE"

# list of markdup BAMs
mapfile -t BAMS < <(find "$BASE_IN" -name "*_realigned_sorted_markdup.bam" | sort)

IN_BAM="${BAMS[$((SLURM_ARRAY_TASK_ID - 1))]}"
SAMPLE="$(basename "$IN_BAM" _realigned_sorted_markdup.bam)"

# Explicit sample identification in logs
echo "SAMPLE=${SAMPLE}"
echo "SAMPLE=${SAMPLE}" >&2

OUT_BED_GZ="${BASE_OUT}/${SAMPLE}.bed.gz"
TMP_DIR="${TMP_BASE}/${SAMPLE}"

mkdir -p "$TMP_DIR"

#No duplicate reads (0x400), secondary alignments (0x100), supplementary alignments (0x800), unmapped reads (0x4)

samtools view \
    -F 3332  \
    -q 20 \
    -b "$IN_BAM" | \
bedtools bamtobed \
    -i - | \
sort -T "$TMP_DIR" -k1,1 -k2,2n | \
bgzip -c > "$OUT_BED_GZ"
```

### Peak calling
Peaks were called using `macs3` version 3.0.1.

```bash
macs2 callpeak \
    -f BED 
    -g hs
    -q 0.05 
    --shift -75 
    --extsize 150 
    --nolambda 
    --nomodel 
    --keep-dup all 
    --call-summits 
```

```bash
#!/bin/bash
#SBATCH --job-name=macs3
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/macs3/%x_%a.out
#SBATCH --error=logs/macs3/%x_%a.err
#SBATCH --array=1-38
#SBATCH --partition=small

set -euo pipefail

# Put temporary files on scratch/job-local storage instead of shared /tmp
if [[ -n "${SLURM_TMPDIR:-}" ]]; then
  export TMPDIR="$SLURM_TMPDIR/macs3_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-0}"
else
  export TMPDIR="../data/tmp/macs3_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-0}"
fi
mkdir -p "$TMPDIR"
echo "TMPDIR=$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT


module load MACS3/3.0.1
macs3 --version

BASE_IN="../data/beds"
BASE_OUT="../data/peaks"

mkdir -p "$BASE_OUT"

# list of BED.gz files
mapfile -t BEDS < <(find "$BASE_IN" -name "*.bed.gz" | sort)
N=${#BEDS[@]}
if (( N == 0 )); then
  echo "No BED.gz files found under $BASE_IN" >&2
  exit 1
fi
if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
  echo "Invalid array task ID ${SLURM_ARRAY_TASK_ID}; valid range is 1..${N}" >&2
  exit 1
fi

IN_BED="${BEDS[$((SLURM_ARRAY_TASK_ID - 1))]}"
SAMPLE="$(basename "$IN_BED" .bed.gz)"

# Sample identification in logs
echo "SAMPLE=${SAMPLE}"
echo "SAMPLE=${SAMPLE}" >&2

macs3 callpeak \
    -t "$IN_BED" \
    --name "$SAMPLE" \
    -f BED \
    -g hs \
    --outdir "$BASE_OUT" \
    -q 0.05 \
    --shift -75 \
    --extsize 150 \
    --nolambda \
    --nomodel \
    --keep-dup all \
    --call-summits
```

### Calling consensus peaks
Consensus peaks across samples were called using `createIterativeOverlapPeakSet.R` ([Repository](https://github.com/corceslab/ATAC_IterativeOverlapPeakMerging)) and the following parameters.

```bash
module load R
Rscript --version # Rscript (R) version 4.4.1 (2024-06-14)

Rscript createIterativeOverlapPeakSet.R \
    --blacklist ENCFF356LFX_and_hg38.blacklist.bed 
    --genome hg38 
    --spm 5 
    --rule "2" 
    --extend 250
```

#### How the exclusion set was downloaded
```bash
# from: https://github.com/Boyle-Lab/Blacklist?tab=readme-ov-file
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz

for f in *.gz; do gunzip $f; done

cat ENCFF356LFX.bed hg38.blacklist.bed | bedtools sort | bedtools merge > ENCFF356LFX_and_hg38.blacklist.bed
```

### Quantification/background correction (part 1) 
Whole-genome extraction of read counts used for background correction. Not counting (i.e., removing) **duplicate reads**, **secondary alignments**, **supplementary alignments**, and **unmapped reads**, as well as reads wth `MAPQ < 20`.

```bash
#!/bin/bash
#SBATCH --job-name=read_counts
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000
#SBATCH --output=logs/read_counts/%x_%a.out
#SBATCH --error=logs/read_counts/%x_%a.err
#SBATCH --array=1-38
#SBATCH --partition=small

set -euo pipefail

module load SAMtools
samtools --version

module load BEDTools
bedtools --version

BASE_IN="../data/bams_markdup"
BASE_OUT="../data/read_counts"

GENOME_FILE="../references/GRCh38/hg38.chrom.sizes.sorted.txt"
GRID_BED="../data/grid_w500_s250.bed.gz"

mkdir -p "$BASE_OUT"

# List of markdup BAMs
mapfile -t BAMS < <(find "$BASE_IN" -name "*_realigned_sorted_markdup.bam" | sort)

IN_BAM="${BAMS[$((SLURM_ARRAY_TASK_ID - 1))]}"
SAMPLE="$(basename "$IN_BAM" _realigned_sorted_markdup.bam)"

# Sample identification in logs
echo "SAMPLE=${SAMPLE}"
echo "SAMPLE=${SAMPLE}" >&2

OUT_TSV_GZ="${BASE_OUT}/${SAMPLE}.tsv.gz"

samtools view \
    -F 3332 \
    -q 20 \
    -b "$IN_BAM" | \
samtools sort -n - | \
bedtools bamtobed \
    -bedpe \
    -i - | \
awk 'BEGIN{OFS="\t"} $1==$4 {print $1,$2,$6}' | \
sort -k1,1 -k2,2n | \
bedtools coverage \
    -counts \
    -sorted \
    -g "$GENOME_FILE" \
    -a "$GRID_BED" \
    -b - | \
gzip > "$OUT_TSV_GZ"
```

### Quantification/background correction (part 2) 
Background correction was performed for the whole-genome.

```bash
#!/bin/bash
#SBATCH --job-name=bg_correct
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/bg_correct/%x_%A_%a.out
#SBATCH --error=logs/bg_correct/%x_%A_%a.err
#SBATCH --array=1-38
#SBATCH --partition=small

set -euo pipefail

# directory with input files
INDIR=../data/read_counts

# Deterministic list of input TSVs
mapfile -t TSVS < <(find "$INDIR" -maxdepth 1 -name "*.tsv.gz" | sort)

# select file for this task
SAMPLE_FILE="${TSVS[$((SLURM_ARRAY_TASK_ID - 1))]}"

# Explicit sample identification in logs
SAMPLE="$(basename "$SAMPLE_FILE" .tsv.gz)"

echo "SAMPLE=${SAMPLE}"
echo "SAMPLE=${SAMPLE}" >&2
echo "Processing: ${SAMPLE_FILE}"

python correct_background.py "${SAMPLE_FILE}" 25
```

```py
import numpy as np
import sys
import glob
import os


def zopen(fn, mode='rt'):
    if fn.endswith('.gz'):
        import gzip
        f = gzip.open(fn, mode)
    else:
        f = open(fn, mode)
    return f


sample_path = sys.argv[1]
sample_name = os.path.basename(sample_path)
percentile = int(sys.argv[2])

gfile_p = '../references/GRCh38/hg38.chrom.sizes.sorted.txt'
centro_p = '../data/centromere_position.tsv'
out_dir = '../data/read_counts_bg_corrected/'
suffix = 'bg_corrected.tsv.gz'  # e.g., bg_corrected.tsv.gz
# out_file = '%s%s_%s' % (out_dir, sample_name, suffix)
out_file = f"{out_dir}{sample_name.replace('.tsv.gz','')}_{suffix}"

step = 250  # CHANGE THIS TO APPROPRAITE VALUES
start = 250  # CHANGE THIS TO APPROPRAITE VALUES

########################
## loading the sample ##
########################

acceptable_chroms = ['chr%s' % i for i in range(1, 23)] + ['chrX', 'chrY']

sample = {}
f = zopen(sample_path, mode='rt')
for l in f:
    chrom, start, end, val = l.rstrip().split('\t')
    if chrom not in acceptable_chroms:
        continue
    sample.setdefault(chrom, [])
    sample[chrom].append(int(val))
f.close()

##################################
## extracting the genome length ##
##################################

gfile = {}
f = open(gfile_p, 'r')
for l in f:
    chrom, length = l.rstrip().split('\t')
    if chrom not in acceptable_chroms:
        continue
    gfile[chrom] = int(length)
f.close()

######################################################################
## extracting the position of the centromeres  key:[p_end, q_start] ##
######################################################################

centro_tmp = {}
f = open(centro_p, 'r')
for l in f:
    chrom, start, end, pq, whereabout = l.rstrip().split('\t')
    if chrom not in acceptable_chroms:
        continue
    centro_tmp.setdefault(chrom, [])
    centro_tmp[chrom].append(int(start))
    centro_tmp[chrom].append(int(end))
f.close()

centro = {}
for chrom in centro_tmp:
    centro[chrom] = [min(centro_tmp[chrom]), max(centro_tmp[chrom])]

del centro_tmp

######################################
## finding the chrom arm bounderies ##
######################################

chrom_arm_boundaries = {}
for chrom in gfile:
    chrom_arm_boundaries[chrom] = {'p': [0, min(centro[chrom])], 'q': [
        max(centro[chrom]), gfile[chrom]]}


################################################
## calculating the chromosome arm percentiles ##
################################################

chrom_arm_percentiles = {}

for chrom, boundaries in chrom_arm_boundaries.items():
    arms = {}
    for arm in boundaries:
        start, end = boundaries[arm]
        arms[arm] = np.percentile(
            sample[chrom][start//step:end//step], percentile)

    chrom_arm_percentiles[chrom] = arms

#################
## corrections ##
#################

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

f = zopen(out_file, 'wt')

## DONOT MOVE THIS LINE!
start, step = 250, 250 ## CHANGE THIS TO APPROPRAITE VALUES

# for chrom in sample:
#     for i, v in enumerate(sample[chrom]):

for chrom, values in sample.items():
    for i, v in enumerate(values):

        current_window_start, current_window_end = (
            i+1) * start - step, (i+1) * start + step

        tk = sample[chrom][0 if (i - 19) < 0 else (i - 19):(i + 2 + 18)]  # 10k
        hk = sample[chrom][0 if (i - 199) < 0 else (i - 199):(i + 2 + 198)]  # 100k

        ###########################################################################
        ## removing the value for the current window and finding the percentiles ##
        ###########################################################################

        tk.remove(v)
        hk.remove(v)

        #######################################
        ## exctracting chrom arm percentiles ##
        #######################################

        pstart, pend = chrom_arm_boundaries[chrom]['p']
        qstart, qend = chrom_arm_boundaries[chrom]['q']

        if current_window_start in range(pstart, pend):  # i.e. window in p arm
            ar = chrom_arm_percentiles[chrom]['p']
        elif current_window_start in range(qstart, qend):
            ar = chrom_arm_percentiles[chrom]['q']
        elif current_window_start < 0:  # this is for the beginning of each chromosome
            ar = chrom_arm_percentiles[chrom]['p']

        bgs = [np.percentile(tk, percentile),
               np.percentile(hk, percentile), ar]
        bgs_str = ', '.join([str(bg) for bg in bgs])

        max_correction_val = max(bgs)
        corrected_val = v - max_correction_val
        corrected_val = 0 if corrected_val < 0 else corrected_val

        f.write('%s\t%s\t%s\t%.2f\t%s\t%.2f\t%s\n' % (chrom, current_window_start,
                current_window_end, corrected_val, v, max_correction_val, bgs_str))

f.close()
```
### Quantification/background correction (part 3) 
A whole-genome matrix of backgrounds was generated

```py
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

%%time 
p = "../data/read_counts_bg_corrected"
files = sorted([f"{p}/{f}" for f in os.listdir(p) if f.endswith('.gz')])
num_files = len(files)

master_df = pd.read_csv(files[0], sep='\t', usecols=[0, 1, 2], header=None)
num_rows = master_df.shape[0]

data_matrix = np.empty((num_rows, num_files), dtype='float32')

for i in range(num_files):
    # Load ONLY column 5 (the data) as float32
    data_matrix[:, i] = pd.read_csv(files[i], sep='\t', header=None, 
                                   usecols=[5], dtype='float32').values.flatten()
    print(f"{i}: Loaded data from {os.path.basename(files[i])}")

master_df.columns = ['chrom', 'start', 'end']

for i in range(num_files):
    sample_name = os.path.basename(files[i]).replace('_bg_corrected.tsv.gz', '')
    master_df[sample_name] = data_matrix[:, i]

## This many take some 10 minutes to write to file
header_str = '\t'.join(master_df.columns)

p = "../data/read_counts_consensus/background.tsv.gz"
np.savetxt(
    p, 
    master_df.values, 
    fmt='%s', 
    delimiter='\t', 
    header=header_str, 
    comments=''
)
```
### Quantification/background correction (part 4) 
A matrix of read counts for the consensus peaks was created using `featureCounts` (v2.1.1).

```bash
#!/bin/bash
#SBATCH --job-name=featureCounts_consensus
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1:30:00
#SBATCH --partition=small
#SBATCH --output=logs/featureCounts/featureCounts_consensus.%j.out
#SBATCH --error=logs/featureCounts/featureCounts_consensus.%j.err

set -euo pipefail

export PATH="/home/USER/local_tools/subread-2.1.1-Linux-x86_64/bin:$PATH"
featureCounts -v

BASE_IN="../data/bams_markdup"
BASE_OUT="../data/read_counts_consensus"
TMP_DIR="${BASE_OUT}/tmp"
SAF_FILE="../data/consensus_peaks/All_Samples.fwp.filter.non_overlapping.SAF.gz"

mkdir -p "${BASE_OUT}"
mkdir -p "${TMP_DIR}"

# List of BAMs files
mapfile -t BAMS < <(find "${BASE_IN}" -name "*_realigned_sorted_markdup.bam" | sort)

echo "Number of BAM files: ${#BAMS[@]}"
echo "BAM files:"
printf '  %s\n' "${BAMS[@]}"

featureCounts \
  -F SAF \
  -a "${SAF_FILE}" \
  -o "${BASE_OUT}/consensus_read_counts.txt" \
  -s 0 \
  -f \
  -O \
  --minOverlap 1 \
  -Q 20 \
  --primary \
  --ignoreDup \
  -p \
  --countReadPairs \
  -T ${SLURM_CPUS_PER_TASK} \
  --tmpDir "${TMP_DIR}" \
  --verbose \
  "${BAMS[@]}"


echo "featureCounts finished: $(date)"
```

#### How the SAF file (used above) was generated

```bash
cd ../data/consensus_peaks

#All_Samples.fwp.filter.non_overlapping.bed is the output of `createIterativeOverlapPeakSet.R`.
awk 'BEGIN{FS="\t"; OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"} {print $1":"$2"-"$3, $1, $2+1, $3, "."}' <(cat All_Samples.fwp.filter.non_overlapping.bed) | gzip > All_Samples.fwp.filter.non_overlapping.SAF.gz
```


- Clean up the `featureCounts` output

```bash
cat ../data/read_counts_consensus/consensus_read_counts.txt | \
awk '
BEGIN {
    OFS = "\t"
    header_seen = 0
}

# Skip all featureCounts metadata lines
/^#/ { next }

# First non-comment line = header
header_seen == 0 {
    header_seen = 1

    # Clean sample names
    for (i = 7; i <= NF; i++) {
        gsub(/^.*\//, "", $i)
        gsub(/_realigned_sorted_markdup\.bam$/, "", $i)
    }

    # Print cleaned header
    printf "chrom\tstart\tend"
    for (i = 7; i <= NF; i++) printf "\t%s", $i
    printf "\n"
    next
}

# Data lines
{
    printf "%s\t%d\t%s", $2, $3 - 1, $4
    for (i = 7; i <= NF; i++) printf "\t%s", $i
    printf "\n"
}
' | gzip > ../data/read_counts_consensus/consensus_read_counts_cleaned.tsv.gz
```

### Quantification/background correction (part 5) 
Consensus peaks were background corrected.  

```py
p = '../data/read_counts_consensus/consensus_read_counts_cleaned.tsv.gz'
consensus_raw = pd.read_csv(p, sep='\t', header=0)
consensus_raw.rename(columns={'chrom':'Chromosome', 'start':'Start', 'end':'End'}, inplace=True)
consensus_raw_pr = pr.PyRanges(consensus_raw)

p = '../data/read_counts_consensus/background.tsv.gz'
backgrounds = pd.read_csv(p, sep='\t', header=0)

backgrounds.rename(columns={'chrom':'Chromosome', 'start':'Start', 'end':'End'}, inplace=True)
backgrounds_pr = pr.PyRanges(backgrounds)

subset = consensus_raw_pr.join(backgrounds_pr)
subset_df = subset.df.drop(columns=['Start_b', 'End_b'])

subset_df = subset_df.groupby(['Chromosome', 'Start', 'End'], observed=True).mean().map(np.rint).astype(int)

cols = [c for c in subset_df.columns if not c.endswith('_b') and f"{c}_b" in subset_df.columns]
for c in cols:
    subset_df[c] = subset_df[c] - subset_df[f"{c}_b"]

subset_df = subset_df.drop(columns=[c for c in subset_df.columns if c.endswith('_b')])

subset_df = subset_df.clip(lower=0)

p = '../data/read_counts_consensus/consensus_read_counts_bg_corrected.tsv.gz'
subset_df.to_csv(p, sep='\t')
``` 
