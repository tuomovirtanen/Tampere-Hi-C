


outpath="../results/corrected_peaks"
peakfile="../results/counts/peak_count_matrix.tsv"
peak_gc_file="../results/peak_gc_fraction.bed"
peak_logratio_file="../results/coverage/peak_copy_logratios.tsv"
coverage_path="../results/coverage"
black_list_file="../results/blacklist_peaks.bed"
reference_pattern="BPH"

#samples to exclude from the analysis due to low coverage etc. samples provided by Ebrahim
exclude_samples=c("PC_15194","PC_20873","PC_22392","PC_22603","PC_24173","PC_4786","PC_6342")

library(edgeR)
library(limma)
library(locfit)
library(Cairo)
library(GenomicRanges)
library(gridExtra)
source("genome_plot.r")
# ---------------------------------------------------------------
# 1. Read and prepare input
# ---------------------------------------------------------------
peaks <- read.table(peakfile,header=TRUE)
colnames(peaks)[1:3] <- c("chrom","start","end")
gc <- read.table("../results/peak_gc_fraction.bed")[,4] #make sure the order is the same as in peakfile


if(length(gc) != nrow(peaks)){
  stop("Error: length of 'gc' does not match number of rows in 'peakfile'")
}


peak_gr <- makeGRangesFromDataFrame(peaks,
                                    keep.extra.columns=FALSE,
                                    ignore.strand=TRUE,
                                    seqnames.field="chrom",
                                    start.field="start",
                                    end.field="end")
seqlevelsStyle(peak_gr)<- "UCSC"
peak_logratios <- read.table(peak_logratio_file,header=TRUE)


black_list <- read.table(black_list_file,header=FALSE)
colnames(black_list)[1:3] <- c("chrom","start","end")
black_list_gr <- makeGRangesFromDataFrame(black_list,
                                    keep.extra.columns=FALSE,
                                    ignore.strand=TRUE,
                                    seqnames.field="chrom",
                                    start.field="start",
                                    end.field="end")
seqlevelsStyle(black_list_gr)<- "UCSC"

keep <- !overlapsAny(peak_gr,black_list_gr)

peaks <- peaks[keep,]
counts <- peaks[,4:ncol(peaks)]
peaks <- peaks[ ,1:3]
gc <- gc[keep]
peak_gr <- peak_gr[keep]
peak_logratios <- peak_logratios[keep,]


samples <- intersect(colnames(peak_logratios), colnames(counts))
samples <- samples[!samples %in% exclude_samples]
peak_logratios <- peak_logratios[, samples]
counts <- counts[, samples]

# ---------------------------------------------------------------
# 2. Set up edgeR object 
# ---------------------------------------------------------------

groups <- factor(gsub("_.*", "", samples)) # or define manually
dge <- DGEList(counts = counts[, samples], group = groups)

#filter the lowly expressed peaks before calculating the offsets
keep <- filterByExpr(dge, group=groups)
peaks <- peaks[keep,]
gc <- gc[keep]
peak_gr <- peak_gr[keep]
peak_logratios <- peak_logratios[keep,]
dge <- dge[keep,]
dge <- calcNormFactors(dge)

design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
rownames(design) <- samples




tss <- read.table("tss_enrcihment.txt",header=TRUE,sep="\t")
all(samples %in% tss$Sample)
all(tss$Sample[match(samples,tss$Sample)] == samples)
tss <- tss$TSS_enrichment[match(samples,tss$Sample)]
design <- cbind(design, tss)
colnames(design)[ncol(design)] <- "tss_enrichment"


# Identify baseline (e.g. BPH) samples
reference <- grepl(reference_pattern, colnames(counts))
reference_baseline <- apply(counts[, reference], 1, median) + 0.1  # avoid division by zero

cp <- cpm(dge, log = TRUE, normalized.lib.sizes=TRUE)
reference_baseline <- apply(cp[,reference],1,median)
log2fold <- apply(cp,2,function(x) x - reference_baseline)

# ---------------------------------------------------------------
# 3. Compute offsets per sample via local regression
# ---------------------------------------------------------------
span <- 0.3
maxk <- 500
offsets <- matrix(NA, nrow = nrow(log2fold), ncol = ncol(log2fold))
rownames(offsets) <- rownames(log2fold)
colnames(offsets) <- colnames(log2fold)

for (i in seq_len(ncol(log2fold))) {
  sample <- colnames(log2fold)[i]
  cat("Fitting offsets for sample:", sample, "\n")
  
  y <- log2fold[, i]                # observed ATAC log2 fold-change
  x <- peak_logratios[, i]          # copy-number log-ratio
  
  # remove missing / infinite values
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  
  # fit local regression (robust) GC bias doesn't have interaction with coverage log2ratio as we have normalized that before when calculating log2fold
  fit <- locfit(y ~ lp(x, nn = span, deg = 1) + lp(gc[keep], nn = span, deg = 1), maxk = maxk, lfproc = locfit.robust)
  
  # predict fitted bias (same order)
  fitted_vals <- rep(NA, length(y))
  fitted_vals[keep] <- fitted(fit)
  
  offsets[, i] <- fitted_vals
}

for(i in 1:ncol(log2fold)){
  s=samples[i]
  print(s)
  #CairoPNG(paste0(outpath,"/offsets_ATAC_copybias_",s,".png"),height=800,width=800)
  pdf(paste0(outpath,"/offsets_ATAC_copybias_",s,".pdf"),height=8,width=8)
    par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))
    plot(peak_logratios[, i], log2fold[, i], pch=20, cex=0.4,
        main=colnames(log2fold)[i], xlab="Copy-number log2 ratio", ylab="ATAC log2 fold-change")
        abline(0, 1, col="black", lty=2) 
    points(peak_logratios[, i], offsets[, i], col="red", pch=20, cex=0.4)
    abline(h=0, col="gray")

    plot(gc, log2fold[, i], pch=20, cex=0.4,
      main=colnames(log2fold)[i], xlab="GC fraction", ylab="ATAC log2 fold-change")
    points(gc, offsets[, i], col="red", pch=20, cex=0.4)
    abline(h=0, col="gray")
  dev.off()
}

# ---------------------------------------------------------------
# 4. Assign offsets and re-estimate dispersion
# ---------------------------------------------------------------

offsets_ln <- offsets * log(2) #Convert to natural log scale (edgeR expects ln offsets)
o <- outer(rep(1, nrow(dge)), getOffset(dge)) + offsets_ln
write.table(o,paste0(outpath,"/offset_matrix.tsv"),sep="\t",quote=FALSE,row.names=TRUE)

#o <- as.matrix(read.table(paste0(outpath,"/offset_matrix.tsv")))
dge$offset <- o

#keep only peaks for which we could estimate the offsets. These were probably problematic regions anyway with low coverage
keep <- rowSums(!is.finite(o)) == 0
dge <- dge[keep,]

dge <- estimateDisp(dge, design = design)
cp <- cpm(dge, log = TRUE, normalized.lib.sizes=TRUE)
colnames(cp) <- samples
write.table(cbind(peaks[keep,],cp),paste0(outpath,"/peaks_cpm_corrected.tsv"),sep="\t",quote=FALSE,row.names=TRUE)

reference_baseline_corrected <- apply(cp[,reference],1,median)
log2fold_corrected <- apply(cp,2,function(x) x - reference_baseline_corrected)


#genome wide plots to estimate the preformance of the normalization
for(s in samples){
  print(s)
  gg <- genome_plot(peaks[keep,1], rowMeans(peaks[keep,2:3]), log2fold[keep,s],ylab="ATAC-log2-ratio") + ylim(quantile(log2fold[keep,s], 0.01), quantile(log2fold[keep,s], 0.99))
  gg2 <- genome_plot(peaks[keep,1], rowMeans(peaks[keep,2:3]), log2fold_corrected[,s],ylab="ATAC-log2-ratio-corrected")
  #png(paste0(outpath,"/",s,".png"),width=10,height=4,units = "in",res=300,type="cairo")
  pdf(paste0(outpath,"/",s,".pdf"),width=10,height=4)
  print(grid.arrange(gg, gg2, ncol = 1))
  dev.off()
}

# ---------------------------------------------------------------
# 5. DE peak analysis
# ---------------------------------------------------------------

#contrast.matrix = eval(as.call(c(as.symbol("makeContrasts"),as.list(myContrasts),levels=list(design))))
#myContrasts <- c("PC-BPH","CRPC-BPH","CRPC-PC")
contrast.matrix <- makeContrasts(
    PC_BPH       = PC - BPH,
    CRPC_BPH     = CRPC - BPH,
    CRPC_PC      = CRPC - PC,
    CRPC_U   = CRPC - (PC + BPH)/2,
    Cancer_BPH   = (PC + CRPC)/2 - BPH,
    levels = design
)


dgefit <- glmQLFit(dge,design=design)
results <- list()
dge_results <- list()
for(i in 1:ncol(contrast.matrix)){
	dge_results[[i]] <- glmQLFTest(dgefit,contrast=contrast.matrix[,i])
	results[[i]] <- dge_results[[i]]$table
	results[[i]]$fdr <- p.adjust(results[[i]]$PValue,method="fdr")
	colnames(results[[i]]) <- paste0(colnames(contrast.matrix)[i],"-",colnames(results[[i]]))
}

results <- do.call(cbind,results)
out <- cbind(peaks[keep,], results)


write.table(out,paste0(outpath,"/diff_peak_results.tsv"),sep="\t",quote=FALSE,row.names=FALSE)


##############

design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
rownames(design) <- samples
tss <- read.table("tss_enrcihment.txt",header=TRUE,sep="\t")
all(samples %in% tss$Sample)
all(tss$Sample[match(samples,tss$Sample)] == samples)
tss <- tss$TSS_enrichment[match(samples,tss$Sample)]
design <- cbind(design, tss)
colnames(design)[ncol(design)] <- "tss_enrichment"

purity_data <- read.table("ATAC-sample-purity.tsv",header=TRUE,sep="\t")
m=match(samples,purity_data$Sample)
tumor_purity <- purity_data$tumor_fraction[m]
tumor_purity[is.na(tumor_purity)] <- 0
tumor_purity
design[!reference,1:3] <- design[!reference,1:3] * (tumor_purity[!reference] / 100)

design[!reference,1] <- (1 - apply(design[!reference,2:3],1,max))
dge <- estimateDisp(dge, design = design)

dgefit <- glmQLFit(dge,design=design)
results <- list()
dge_results <- list()
for(i in 1:ncol(contrast.matrix)){
	dge_results[[i]] <- glmQLFTest(dgefit,contrast=contrast.matrix[,i])
	results[[i]] <- dge_results[[i]]$table
	results[[i]]$fdr <- p.adjust(results[[i]]$PValue,method="fdr")
	colnames(results[[i]]) <- paste0(colnames(contrast.matrix)[i],"-",colnames(results[[i]]))
}
results <- do.call(cbind,results)
out <- cbind(peaks[keep,], results)

write.table(out,paste0(outpath,"/diff_peak_results_tumor_purity.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
#############


design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
rownames(design) <- samples

purity_data <- read.table("ATAC-sample-purity.tsv",header=TRUE,sep="\t")
m=match(samples,purity_data$Sample)
tumor_purity <- purity_data$tumor_fraction[m]
tumor_purity[is.na(tumor_purity)] <- 0
tumor_purity
design[!reference,1:3] <- design[!reference,1:3] * (tumor_purity[!reference] / 100)

design[!reference,1] <- (1 - apply(design[!reference,2:3],1,max))
dge <- estimateDisp(dge, design = design)

dgefit <- glmQLFit(dge,design=design)
results <- list()
dge_results <- list()
for(i in 1:ncol(contrast.matrix)){
	dge_results[[i]] <- glmQLFTest(dgefit,contrast=contrast.matrix[,i])
	results[[i]] <- dge_results[[i]]$table
	results[[i]]$fdr <- p.adjust(results[[i]]$PValue,method="fdr")
	colnames(results[[i]]) <- paste0(colnames(contrast.matrix)[i],"-",colnames(results[[i]]))
}
results <- do.call(cbind,results)
out <- cbind(peaks[keep,], results)

write.table(out,paste0(outpath,"/diff_peak_results_tumor_purity-no-tss.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
#############
