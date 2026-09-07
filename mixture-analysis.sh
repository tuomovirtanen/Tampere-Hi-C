#!/bin/bash -l
#SBATCH -J mix
#SBATCH -o mix_%A_%a.out
#SBATCH -e mix_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH -t 48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -p large
#SBATCH --array=1-10

module load Java

TMP="mixture_analysis/tmp"
mkdir -p "$TMP"

fragment_file1="../fragment_files/BPH456_merged_nodups.txt.gz"
fragment_file2="../fragment_files/CRPC278_merged_nodups.txt.gz"

N1=1278174832
N2=1539257664
TARGET=500000000

FILE1=$fragment_file1
FILE2=$fragment_file2

frac2_values=(100 90 80 70 60 50 40 30 20 10)
FRAC2=${frac2_values[$SLURM_ARRAY_TASK_ID-1]}
FRAC1=$((100 - FRAC2))

P1=$(awk "BEGIN{printf \"%.12f\", ($FRAC1/100.0)*$TARGET/$N1}")
P2=$(awk "BEGIN{printf \"%.12f\", ($FRAC2/100.0)*$TARGET/$N2}")

SORTED="mix_${FRAC2}.merged_nodups.short.sorted.txt.gz"
HIC="mix_${FRAC2}.hic"

(
    zcat "$FILE1" |
    awk -v p="$P1" 'BEGIN{srand()}
        rand()<p && $3==$7 {
            print $2, $3, $4, $5, $6, $7, $8, $9
        }'

    zcat "$FILE2" |
    awk -v p="$P2" 'BEGIN{srand()}
        rand()<p && $3==$7 {
            print $2, $3, $4, $5, $6, $7, $8, $9
        }'
) |
sort -T "$TMP" --parallel=$SLURM_CPUS_PER_TASK \
     -k2,2 -k6,6 -k3,3n -k7,7n |
gzip > "$SORTED"

java -Xms64G -Xmx64G \
    -jar ../juicer_tools_1.14.08.jar \
    pre \
    -d \
    -r 100000 \
    "$SORTED" \
    "$HIC" \
    hg38



###############################
#!/bin/bash -l
#SBATCH -J mix
#SBATCH -o mix_%A_%a.out
#SBATCH -e mix_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH -t 48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -p large
#SBATCH --array=1-10

module load Java

fragment_file1="../fragment_files/BPH456_merged_nodups.txt.gz"
fragment_file2="../fragment_files/CRPC278_merged_nodups.txt.gz"

N1=1278174832
N2=1539257664
TARGET=1000000000
# Set these before submitting or export them in your environment
FILE1=${fragment_file1}
FILE2=${fragment_file2}

TARGET=${TARGET}
N1=${N1}
N2=${N2}

# FRAC2 values corresponding to array IDs 1-10
frac2_values=(100 90 80 70 60 50 40 30 20 10)
FRAC2=${frac2_values[$SLURM_ARRAY_TASK_ID-1]}
FRAC1=$((100 - FRAC2))

P1=$(awk "BEGIN{printf \"%.12f\", ($FRAC1/100.0)*$TARGET/$N1}")
P2=$(awk "BEGIN{printf \"%.12f\", ($FRAC2/100.0)*$TARGET/$N2}")
SORTED=mix_${FRAC2}.merged_nodups.sorted.txt.gz
HIC=mix_${FRAC2}.hic

(
    zcat "$FILE1" |
    awk -v p="$P1" 'BEGIN{srand()} rand()<p && $3==$7'

    zcat "$FILE2" |
    awk -v p="$P2" 'BEGIN{srand()} rand()<p && $3==$7'
) | sort -k3,3 -k7,7 -k4,4n -k8,8n \
  | gzip > "$SORTED"

java -Xms64G -Xmx64G \
    -jar juicer_tools_1.14.08.jar \
    pre \
    -d \
    -r 100000 \
    "$SORTED" \
    "$HIC" \
    hg38


awk 'NR>1 {gsub(/"/,""); print $3"\t"$5"\t"$6}' CRPC278_segments.txt > segments.bed
launch_R
Rscript normalize-hic-segments.r sbatch/mix_100.hic segments.bed test.txt


####generate pca files from the .hic files

#!/bin/bash -l
#SBATCH -J pca
#SBATCH -o pca_%A_%a.out
#SBATCH -e pca_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH -p large
#SBATCH --mem=32G
#SBATCH --array=1-10

set -euo pipefail

samples=(
    mix_100
    mix_10
    mix_20
    mix_30
    mix_40
    mix_50
    mix_60
    mix_70
    mix_80
    mix_90
)

sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

module load Java

cd mixture_analysis/normed

# Remove rows containing NA
grep -v 'NA' "${sample}.txt" > "${sample}_noNA.txt"

# Generate .hic
java -Xmx32G -jar ../juicer_tools_1.14.08.jar \
    pre "${sample}_noNA.txt" "${sample}.hic" hg38 -r 100000 -d

# Compute PC1 for each chromosome
chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

mkdir -p ../pca

for chr in "${chrs[@]}"; do
    echo "Processing ${sample}, chr${chr}"

    java -Xmx32G -jar ../juicer_tools_1.14.08.jar \
        eigenvector -p KR \
        "${sample}.hic" \
        "${chr}" \
        BP \
        100000 \
        "../pca/${sample}-pc1-${chr}.txt"
done


########################

library(GenomicRanges)

samples <- c(
  "mix_100",
  "mix_10",
  "mix_20",
  "mix_30",
  "mix_40",
  "mix_50",
  "mix_60",
  "mix_70",
  "mix_80",
  "mix_90"
)

chrs <- c(as.character(1:22), "X", "Y")
res <- 100000

path <- "mixture_analysis/pca"
outpath <- "mixture_analysis/"


RNA_quant <- read.table(
    "Hi-C/GE_deseq_norm.txt",
    header = TRUE
)

rna_gr <- makeGRangesFromDataFrame(
    RNA_quant,
    keep.extra.columns = FALSE,
    ignore.strand = TRUE,
    seqnames.field = "chromosome_name",
    start.field = "tss",
    end.field = "tss"
)

flip <- function(a, b) {
    ifelse(a >= b, 1, -1)
}

PCs <- list()
bins <- list()

for (chr in chrs) {

    message("Processing chr", chr)

    pc <- lapply(samples, function(s) {

        f <- file.path(path, paste0(s, "-pc1-", chr, ".txt"))

        if (!file.exists(f))
            stop("Missing file: ", f)

        x <- scan(f, quiet = TRUE)

        if (length(x) == 0)
            stop("Empty file: ", f)

        x
    })

    lim <- min(sapply(pc, length))

    pc <- lapply(pc, function(x) x[seq_len(lim)])

    b <- data.frame(
        chr = chr,
        start = (seq_len(lim) - 1) * res,
        end = seq_len(lim) * res
    )

    b_gr <- makeGRangesFromDataFrame(
        b,
        keep.extra.columns = FALSE,
        ignore.strand = TRUE,
        seqnames.field = "chr",
        start.field = "start",
        end.field = "end"
    )

    pc <- lapply(pc, function(x) {

        pos <- x > 0 & !is.na(x)
        neg <- x < 0 & !is.na(x)

        pos.expr <- sum(rna_gr %over% b_gr[pos]) / sum(pos)
        neg.expr <- sum(rna_gr %over% b_gr[neg]) / sum(neg)

        x * flip(pos.expr, neg.expr)
    })

    PCs[[chr]] <- do.call(cbind, pc)
    bins[[chr]] <- b
}

PC <- do.call(rbind, PCs) * 100
colnames(PC) <- samples

Bins <- do.call(rbind, bins)
options(scipen = 999)
pc_out <- cbind(Bins, PC)

write.table(
    pc_out,
    file.path(outpath, "PC1-downsampled.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

########

PC1 <- read.table("Hi-C/PC1.txt",header=TRUE)

PC1_downsampled <- read.table(
    "mixture_analysis/PC1-downsampled.txt",
    header = TRUE
)

all(PC1[,2] == PC1_downsampled[,2])


bins = PC1[,1:3]
PC1 = PC1[,-c(1:3)]


keep=rowSums(is.na(PC1)) == 0 & bins[,1] != "X" & bins[,1] != "Y" & rowSums(is.na(PC1_downsampled)) == 0 
bins<-bins[keep,]
PC1 <- PC1[keep,]
PC1_downsampled <- PC1_downsampled[keep,4:ncol(PC1_downsampled)]


BPH_ref=PC1$BPH456

#save plots to same pdf use ggplot, also fit a line

library(ggplot2)

## Calculate delta relative to BPH mean
delta <- sweep(PC1_downsampled, 1, BPH_ref, "-")

## Samples to compare
samples <- c(
    "mix_100",
    "mix_10",
    "mix_20",
    "mix_30",
    "mix_40",
    "mix_50",
    "mix_60",
    "mix_70",
    "mix_80",
    "mix_90"
)

## Check that all samples exist
stopifnot(all(samples %in% colnames(delta)))
library(MASS)
fit_list <- list()
cairo_pdf("PC1_downsampled_vs_mix100.pdf", width = 5, height = 5)

fractions <- c(
    mix_100 = 1.00,
    mix_10  = 0.10,
    mix_20  = 0.20,
    mix_30  = 0.30,
    mix_40  = 0.40,
    mix_50  = 0.50,
    mix_60  = 0.60,
    mix_70  = 0.70,
    mix_80  = 0.80,
    mix_90  = 0.90
)

for (s in names(fractions)[-1]) {

    df <- data.frame(
        x = delta$mix_100,
        y = delta[[s]]
    )

    fit.robust <- rlm(y ~ x, data = df)
    fit <- lm(y ~ x, data = df)
    fit_list[[s]] <- fit.robust

    r2 <- summary(fit)$r.squared
    corval <- cor(df$x, df$y)

    newx <- data.frame(x = seq(min(df$x), max(df$x), length.out = 200))

    p <- ggplot(df, aes(x, y)) +
        geom_point(alpha = 0.3, size = 0.4) +

        ## ordinary least squares
        geom_abline(
            intercept = coef(fit)[1],
            slope = coef(fit)[2],
            colour = "red",
            linewidth = 1
        ) +

        ## robust regression
        geom_line(
            data = data.frame(
                x = newx$x,
                y = predict(fit.robust, newdata = newx)
            ),
            aes(x, y),
            inherit.aes = FALSE,
            colour = "darkgreen",
            linewidth = 1
        ) +

        ## theoretical dilution line
        geom_abline(
            intercept = 0,
            slope = fractions[s],
            colour = "blue",
            linetype = 2,
            linewidth = 1
        ) +

        labs(
            title = paste("mix_100 vs", s),
            x = expression(Delta*"PC1 (100%)"),
            y = paste0("\u0394PC1 (", s, ")"),
            subtitle = sprintf(
                "r = %.3f   R² = %.3f\nLM slope = %.3f   Robust slope = %.3f   Expected = %.2f",
                corval,
                r2,
                coef(fit)[2],
                coef(fit.robust)[2],
                fractions[s]
            )
        ) +
        theme_bw()

        print(p)
    }

dev.off()


cols <- setNames(
    colorRampPalette(c("dodgerblue3", "orange", "firebrick"))(9),
    names(fractions)[-1]
)

xrange <- range(delta$mix_100, na.rm = TRUE)
xseq <- seq(xrange[1], xrange[2], length.out = 200)
cairo_pdf("PC1_downsampled_summary.pdf", width = 5, height = 5)
p <- ggplot() +
    theme_bw() +
    labs(
        title = "PC1 dilution model",
        x = expression(Delta*"PC1 (100%)"),
        y = expression(Delta*"PC1 (mixture)")
    )

for (s in names(fractions)[-1]) {

    ## theoretical line
    p <- p +
        geom_abline(
            intercept = 0,
            slope = fractions[s],
            colour = cols[s],
            linetype = 2,
            linewidth = 0.8
        )

    ## robust fit
    pred <- data.frame(
        x = xseq,
        y = predict(fit_list[[s]], newdata = data.frame(x = xseq))
    )

    p <- p +
        geom_line(
            data = pred,
            aes(x, y),
            colour = cols[s],
            linewidth = 1.2
        )
}

print(p)

dev.off()



library(dplyr)

cols <- setNames(
    colorRampPalette(c("dodgerblue3", "orange", "firebrick"))(9),
    names(fractions)[-1]
)

xrange <- range(delta$mix_100, na.rm = TRUE)
xseq <- seq(xrange[1], xrange[2], length.out = 200)

## Expected lines
expected_df <- do.call(rbind, lapply(names(fractions)[-1], function(s) {
    data.frame(
        x = xseq,
        y = fractions[s] * xseq,
        Fraction = s,
        Type = "Expected"
    )
}))

## Robust fits
robust_df <- do.call(rbind, lapply(names(fractions)[-1], function(s) {
    data.frame(
        x = xseq,
        y = predict(fit_list[[s]], newdata = data.frame(x = xseq)),
        Fraction = s,
        Type = "Robust"
    )
}))

line_df <- rbind(expected_df, robust_df)

pdf("PC1_summary_fits.pdf", width = 6, height = 5)

ggplot(line_df,
       aes(x, y,
           colour = Fraction,
           linetype = Type,
           group = interaction(Fraction, Type))) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = cols) +
    theme_bw() +
    labs(
        title = "PC1 dilution model",
        x = expression(Delta*"PC1 (100%)"),
        y = expression(Delta*"PC1 (mixture)"),
        colour = "Mixture",
        linetype = "Line"
    )+ylim(-5,5)+xlim(-5,5)

dev.off()