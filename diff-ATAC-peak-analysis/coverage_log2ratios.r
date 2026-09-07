count_path="../results/counts"
bins_file="../results/bins_100000.bed"
outpath="../results/coverage"
referecen_pattern="BPH" # pattern to identify reference samples grepl(referecen_pattern,colnames(count_matrix))
gc_file="../results/bins_100000_GC_fraction.bed"
genomeplot="./genome_plot.r"
peakfile="../results/peaks.bed"


library(GenomicRanges)
library(locfit)
library(ggplot2)
library(DNAcopy)

source(genomeplot)

files=list.files(count_path,pattern = "*_bin_counts.bed")
bins=read.table(bins_file)

count_matrix <- list()
for(f in files){
  counts <- read.table(paste0(count_path,"/",f))
  #extract sample name
  s <- sub("_bin_counts.bed.gz$", "", f)
  #pool the masked bins to their original bins
  pooled <- split(counts$V5, counts$V4) #5th bin beacause the bed includes the bin ids in 4th column
  pooled <- unlist(lapply(pooled,sum))
  idx=as.numeric(names(pooled))
  out <- rep(NA,nrow(bins))
  out[idx] <- pooled
  count_matrix[[s]] <- out
}
count_matrix <- do.call(cbind,count_matrix)
out <- cbind(bins,count_matrix)
colnames(out)[1:3] <- c("chr","start","end")
write.table(out,paste0(outpath,"/count_matrix.tsv"),quote=F,sep="\t",col.names=T,row.names=F)

gc=read.table(gc_file)[,4]
#this works well with 100kb bins
extreme_gc <- gc < 0.3 | gc > 0.7

#require that all samples have at least 10 counts in a bin
#this ensures that there are no zero counts when calculating ratios
#if there are regions with deep deletions in high tumor purity samples etc that might be filtered out here.
#if you modify this so that zero counts are allowed e.g. in one sample, make sure to add a pseudocount when calculating ratios
keep <- rowSums(count_matrix > 10) == ncol(count_matrix) & !extreme_gc

#divide each sample by their median
log2ratios <- apply(count_matrix[keep,],2,function(x) x / median(x))

reference <- grepl(referecen_pattern,colnames(count_matrix))
reference_median <- apply(log2ratios[,reference],1,median)

#calculate log2 ratios
log2ratios <- log2(apply(log2ratios,2, function(x) x / reference_median))

gc_corrected <- list()
for(s in colnames(log2ratios)){
  print(s)
  y=log2ratios[,s]
  x=gc[keep]
  fit <- locfit(y ~ lp(x, nn = 0.3, deg = 1), maxk = 500, lfproc = locfit.robust)
  df <- data.frame(x = x, y = y)

  # Generate fitted values on a sorted x-grid
  x_grid <- seq(min(x), max(x), length.out = 500)
  fit_pred <- predict(fit, newdata = data.frame(x = x_grid))
  fit_df <- data.frame(x = x_grid, y = fit_pred)

  gg <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_line(data = fit_df, aes(x = x, y = y), color = "red", linewidth = 1) +  # fitted line
  labs(x = "GC Content", y = "log2 Ratio", 
       title = "Locfit GC Correction") +
  theme_bw()

  pdf(paste0(outpath,"/GC-bias-",s,".pdf"),width=6,height=6)
  print(gg)
  dev.off()

  offset <- predict(fit, newdata = data.frame(x = x))
  gc_corrected[[s]] <- y - offset
}

gc_corrected <- do.call(cbind,gc_corrected)
out <- cbind(bins[keep,1:3],gc_corrected)
colnames(out)[1:3] <- c("chr","start","end")
write.table(out,paste0(outpath,"/logratio_matrix.tsv"),quote=F,sep="\t",col.names=T,row.names=F)

chr=bins[,1][keep]
width=bins[1,3]-bins[1,2] +1
pos=bins[keep,2] + width


segment_list <- list()
for(s in colnames(gc_corrected)){
    CNA.object <- CNA(
      genomdat = gc_corrected[,s],
      chrom = chr,
      maploc = pos,   # or midpoints if you prefer
      data.type = "logratio",
      sampleid = s
    )
  #smoothed.CNA.object <- smooth.CNA(CNA.object, smooth.region = 0,outlier.SD.scale = 4)
  segment.CNA.object <- segment(CNA.object, verbose = 1,  alpha = 0.01,undo.splits = "sdundo",undo.SD = 4,min.width = 2)
  seg_results <- segment.CNA.object$output

  segments=data.frame(
    "chr"=seg_results$chrom,
    "start"= (seg_results$loc.start) - width,
    "end"= seg_results$loc.end + width,
    "log2ratio"= seg_results$seg.mean)

  gg <- genome_plot(chr, pos, gc_corrected[,s],segments=segments,ylab="Coverage log2 ratio")

  write.table(segments,paste0(outpath,"/segmented-logratios-",s,".bed"),quote=F,sep="\t",col.names=T,row.names=F)
  pdf(paste0(outpath,"/logratio-plot-",s,".pdf"),width=10,height=4)
  print(gg)
  dev.off()
  segment_list[[s]] <- segments
}

###logratios for peaks

# ---------------------------------------------------------------
# 1. Read and prepare input
# ---------------------------------------------------------------
peaks <- read.table(peakfile,header=TRUE)
colnames(peaks)[1:3] <- c("chr","start","end")
peak_gr <- makeGRangesFromDataFrame(peaks,
                                    keep.extra.columns=FALSE,
                                    ignore.strand=TRUE,
                                    seqnames.field="chr",
                                    start.field="start",
                                    end.field="end")
seqlevelsStyle(peak_gr)<- "UCSC"


peak_logratios=list()
for(s in names(segment_list)){
  print(s)
  
  cnv_gr <- makeGRangesFromDataFrame(segment_list[[s]],
                                    keep.extra.columns=TRUE,
                                    ignore.strand=TRUE,
                                    seqnames.field="chr",
                                    start.field="start",
                                    end.field="end")
  seqlevelsStyle(cnv_gr)<- "UCSC"
  ols <- findOverlaps(peak_gr,cnv_gr)
  logratios <- rep(NA,length(peak_gr))
  logratios[queryHits(ols)] <- cnv_gr$log2ratio[subjectHits(ols)]
  peak_logratios[[s]] <- logratios
}
peak_logratios = do.call(cbind,peak_logratios)

write.table(peak_logratios,paste0(outpath,"/peak_copy_logratios.tsv"),sep="\t",quote=FALSE,row.names=FALSE)












############
