
library(GenomicRanges)
library(reshape2)

samples <- c("BPH456",  "BPH651",  "BPH659",  "BPH671",  "BPH688",  "BPH701",
  "CRPC278",  "CRPC305",  "CRPC435",  "CRPC489",  "CRPC539",  "CRPC697",
    "PC15420",  "PC15760",  "PC18307",  "PC19403",  "PC4980",  "PC6488")

chrs <- c(1:22,"X","Y")
res=100000
path="."
outpath="."

RNA_quant <- read.table("GE_deseq_norm.txt")
rna_gr <- makeGRangesFromDataFrame(RNA_quant,
                         keep.extra.columns=FALSE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field="chromosome_name",
                         start.field="tss",
                         end.field="tss")



flip <- function(a,b){
    ifelse(a>=b,1,-1)
}

PCs <- list()
bins <- list()

for(chr in chrs){
  pc <- lapply(samples,function(s) read.table(paste0(path,"/",s,"-pc1-",chr,".txt"))[,1])
  lim <- min(unlist(lapply(pc,length)))
  pc <- lapply(pc,function(x) x[1:lim])
  b <- data.frame(chr=chr,start=((1:lim)-1)*res,end=(1:lim)*res)

  b_gr <- makeGRangesFromDataFrame(b,
                         keep.extra.columns=FALSE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field="chr",
                         start.field="start",
                         end.field="end")

  pc <- lapply(pc, function(x) x*flip(sum(rna_gr %over% b_gr[x > 0 & !is.na(x)])/sum(x > 0,na.rm=T),
                                      sum(rna_gr %over% b_gr[x < 0& !is.na(x)])/sum(x < 0,na.rm=T)))

  pc <- do.call(cbind,pc)
  PCs[[chr]] <- pc
  bins[[chr]] <- b
}


PC <- do.call(rbind,PCs)*100
colnames(PC) <- samples
Bins <- do.call(rbind,bins)

pc_out <- cbind(Bins,PC)
write.table(pc_out,paste0(outpath,"/PC1.txt"))

