
#e.g.
#Rscript my_script.R peak_file=/path/to/peaks outpath=/home/user/results

#change these for desired parameters
args <- commandArgs(trailingOnly = TRUE)

# Default values
bin_size <- 100000
padding  <- 500
peak_file <- NULL
outpath <- NULL

# Parse key=value arguments
for (arg in args) {
  key_val <- strsplit(arg, "=")[[1]]
  if (length(key_val) == 2) {
    key <- key_val[1]
    val <- key_val[2]
    if (key == "peak_file") peak_file <- val
    if (key == "outpath")   outpath   <- val
    if (key == "bin_size")  bin_size  <- as.numeric(val)
    if (key == "padding")   padding   <- as.numeric(val)
  }
}

# ---- Check mandatory arguments ----
if (is.null(peak_file) || is.null(outpath)) {
  cat("Error: Missing required arguments.\n\n")
  cat("Usage:\n")
  cat("  Rscript my_script.R peak_file=/path/to/file.tsv.gz outpath=/path/to/output [bin_size=100000] [padding=250]\n\n")
  quit(status = 1)
}

# ---- Print parameters ----
cat("Running with parameters:\n")
cat("  peak_file:", peak_file, "\n")
cat("  outpath  :", outpath, "\n")
cat("  bin_size :", bin_size, "\n")
cat("  padding  :", padding, "\n\n")



library(GenomicRanges)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
options(scipen=999)#don't use scientific notation
#generate xxxkb bins for the whole genome from these chromosome sizes
chrom_lengths <- seqlengths(Seqinfo(genome = "hg38"))
chrom_lengths <- chrom_lengths[grepl("^chr[0-9XY]+$", names(chrom_lengths))]

#generate xxxkb  bins
bins_list <- lapply(names(chrom_lengths), function(chr) {
  chr_len <- chrom_lengths[[chr]]
  starts <- seq(1, chr_len, by = bin_size)
  ends <- pmin(starts + bin_size - 1, chr_len)
  GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends))
})
bins_gr <- do.call(c, bins_list)
#makesure that the seqlevels style (i.e. chromosomes) is UCSC
seqlevelsStyle(bins_gr) <- "UCSC"

#####
library(GenomicRanges)
peaks <- read.table(peak_file,header=FALSE,sep="\t")

#write.table(peaks[,1:3],paste0(outpath,"/peaks.bed"),quote=F,sep="\t",col.names=FALSE,row.names=F)
colnames(peaks)[1:3] <- c("chrom","start","end") #make sure the first 3 columns are chrom,start,end
peak_gr <- makeGRangesFromDataFrame(peaks,
                                    keep.extra.columns=FALSE,
                                    ignore.strand=TRUE,
                                    seqnames.field="chrom",
                                    start.field="start",
                                    end.field="end")

start(peak_gr) <- start(peak_gr) - padding
end(peak_gr) <- end(peak_gr) + padding


#after adding the padding some peaks might overlap, so we need to merge them
peak_gr <- reduce(peak_gr)
#makesure that the seqlevels style (i.e. chromosomes) is UCSC
seqlevelsStyle(peak_gr)<- "UCSC"


new_bins <- list()
for(chr in unique(seqnames(bins_gr))){
  print(chr)
  breaks1 <- c(start(bins_gr[seqnames(bins_gr) == chr]), end(bins_gr[seqnames(bins_gr) == chr]))
  breaks2 <- c(start(peak_gr[seqnames(peak_gr) == chr]), end(peak_gr[seqnames(peak_gr) == chr]))
  breaks <- unique(sort(c(breaks1,breaks2)))
  starts=breaks[1:length(breaks)-1]
  ends=breaks[2:length(breaks)]

  width <- ends - starts
  #we remove bins that are of width 1 because we are just sorting all the breaks and the end of the previous bin and the start of the next bin are always seprated by 1 base
  #this will also remove any potential ultra short bin
  keep <- width > 1
  out <- data.frame("chr"=chr,"start"=starts[keep],ends=ends[keep])
  new_bins[[chr]] <- out
}
new_bins <- do.call(rbind,new_bins)

new_bins_gr <- makeGRangesFromDataFrame(new_bins,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=TRUE,
                                    seqnames.field="chr",
                                    start.field="start",
                                    end.field="ends")

#we now subtract the original peaks (with padding) from these new bins
#but first we need to shrink the peaks by 1 base from each side to avoid removing bins that just touch the peak borders
start(peak_gr) <- start(peak_gr) + 1
end(peak_gr) <- end(peak_gr) - 1

new_bins_gr <- new_bins_gr[!(new_bins_gr %over% peak_gr)]

bin_id <- findOverlaps(new_bins_gr,bins_gr)
bin_id <- subjectHits(bin_id)
out <- data.frame(new_bins_gr)[,1:3]
out$bin_id <- bin_id
write.table(out,paste0(outpath,"/bins_",bin_size,"_masked.bed"),quote=F,sep="\t",col.names=F,row.names=F)
write.table(data.frame(bins_gr)[,1:3],paste0(outpath,"/bins_",bin_size,".bed"),quote=F,sep="\t",col.names=F,row.names=F)


#Calculate GC fraction for each bin
genome <- BSgenome.Hsapiens.UCSC.hg38
# Prepare result vector
gc_fraction <- numeric(length(new_bins_gr))
#I feel this could be faster, but we only need to run this once so it doesn't really matter
# Process chromosome by chromosome
message("Calculating GC fractions for bins...")
for (chr in seqlevels(new_bins_gr)) {
  chr_idx <- which(seqnames(new_bins_gr) == chr)
  if (length(chr_idx) == 0) next
  
  message("Processing ", chr, " (", length(chr_idx), " bins)")
  
  chr_bins <- new_bins_gr[chr_idx]
  chr_seq <- genome[[chr]]  # this is a DNAString object for the full chromosome

  # Create views on the chromosome sequence
  chr_views <- Views(chr_seq, ranges(chr_bins))

  # Compute GC fraction per bin
  gc_fraction[chr_idx] <- rowSums(letterFrequency(chr_views, c("G", "C"), as.prob = TRUE))
}
out_gc <- cbind(out[,1:3],gc_fraction)

write.table(out_gc,paste0(outpath,"/bins_",bin_size,"_masked_GC_fraction.bed"),quote=F,sep="\t",col.names=F,row.names=F)


width <- out_gc$end - out_gc$start + 1
gc_pooled <- split(gc_fraction, bin_id)
width_pooled <- split(width, bin_id)

#calculate GC fraction weighted by width of each masked "sub-bin"
gc_width_pooled <- unlist(lapply(names(width_pooled), function(x) sum(gc_pooled[[x]]*width_pooled[[x]]) / sum(width_pooled[[x]])))

out_gc_bins <- data.frame(data.frame(bins_gr)[,1:3],gc_fraction=gc_width_pooled)
write.table(out_gc_bins,paste0(outpath,"/bins_",bin_size,"_GC_fraction.bed"),quote=F,sep="\t",col.names=F,row.names=F)


#######
#peaks without any padding
peak_original_gr <- makeGRangesFromDataFrame(peaks,
                                    keep.extra.columns=FALSE,
                                    ignore.strand=TRUE,
                                    seqnames.field="chrom",
                                    start.field="start",
                                    end.field="end")


peak_gc_fraction <- numeric(length(peak_original_gr))

# Process chromosome by chromosome
message("Calculating GC fractions for peaks...")
for (chr in seqlevels(peak_original_gr)) {
  chr_idx <- which(seqnames(peak_original_gr) == chr)
  if (length(chr_idx) == 0) next
  chr_bins <- peak_original_gr[chr_idx]
  chr_seq <- genome[[chr]]  # this is a DNAString object for the full chromosome

  # Create views on the chromosome sequence
  chr_views <- Views(chr_seq, ranges(chr_bins))

  # Compute GC fraction per bin
  peak_gc_fraction[chr_idx] <- rowSums(letterFrequency(chr_views, c("G", "C"), as.prob = TRUE))
}

out_gc_peaks <- cbind(peaks[,1:3],peak_gc_fraction)
write.table(out_gc_peaks,paste0(outpath,"/peak_gc_fraction.bed"),quote=F,sep="\t",col.names=F,row.names=F)
