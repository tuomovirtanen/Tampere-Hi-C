
library(GenomicRanges)

PC1 <- read.table("PC1.txt",header=TRUE)
bins = PC1[,1:3]
PC1 = PC1[,-c(1:3)]

keep=rowSums(is.na(PC1)) == 0
bins<-bins[keep,]
PC1 <- PC1[keep,]

bins$chr=paste0("chr",bins$chr)
bins$start <- bins$start +1
bins_gr <- makeGRangesFromDataFrame(bins,
                                   keep.extra.columns=FALSE,
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field="chr",
                                   start.field="start",
                                   end.field="end")
seqlevelsStyle(bins_gr) <- "UCSC"
#define lamin anchor regions
lamin_files <- list.files("data/Lamin-B/chromHMM")
lamin_d <- lapply(lamin_files,function(x) read.table(paste0("data/Lamin-B/chromHMM/",x)))
lamin_d <- lapply(lamin_d,function(x) x[grepl("T1",x[,4]) | grepl("T1",x[,4]),1:3])

lamin_gr <- lapply(lamin_d, function(x) makeGRangesFromDataFrame(x,
                                   keep.extra.columns=FALSE,
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field="V1",
                                   start.field="V2",
                                   end.field="V3"))
ols <- do.call(cbind,lapply(lamin_gr, function(x) bins_gr %over% x))
table(rowSums(ols))

lamin_b_comp <- rowSums(ols) > 9

write.table(data.frame(bins[lamin_b_comp,]),quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE,file="Lamin-anchors.bed")
PC1$BPH_mean <- rowMeans(PC1[,1:6])