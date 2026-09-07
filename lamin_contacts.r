PC1 <- read.table("Hi-C/PC1.txt")
bins = PC1[,1:3]
bins$chr <- paste0("chr",bins$chr)
PC1 = PC1[,-c(1:3)]

lamin <- read.table("Hi-C/Lamin-anchors.bed",header=TRUE)


library(GenomicRanges)
lamin_gr <- makeGRangesFromDataFrame(lamin,
                                   keep.extra.columns=FALSE,
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field="chr",
                                   start.field="start",
                                   end.field="end")
seqlevelsStyle(lamin_gr) <- "UCSC"


hic_file="Hi-C/hic_files/BPH456.hic"
chroms <- c(paste0("chr", 1:22))
chroms
library(Matrix)
library(ggplot2)
for(i in 1:(length(chroms)-1)){
    chr=chroms[i]
    mat <- strawr::straw("SCALE",hic_file,gsub("chr","",chr),gsub("chr","",chr),"BP",100000,"Observed")
    mat$x <- floor(mat$x / 100000) +1
    mat$y <- floor(mat$y / 100000) +1

    mat$counts[is.na(mat$counts)]  <- 0

    dim=min(max(mat$x), max(mat$y))
    mat <- mat[mat$x <= dim & mat$y <= dim,]
    hic_sparse <- sparseMatrix(
      i = mat$x,
      j = mat$y,
      x = mat$counts,
      dims = c(dim,dim),
      symmetric = TRUE # or TRUE for cis-chromosomal maps if appropriate
    )

    chr_bins <- data.frame(chr=chr,start=(seq(1,nrow(hic_sparse)) - 1) * 100000 + 50000,end=(seq(1,nrow(hic_sparse))) * 100000 - 50000)
    chr_bins_gr <- makeGRangesFromDataFrame(chr_bins,
                                   keep.extra.columns=FALSE,
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field="chr",
                                   start.field="start",
                                   end.field="end")


    lamin_idx = chr_bins_gr %over% lamin_gr
    
    lamin_contacts=log10(colSums(hic_sparse[lamin_idx,]) / colSums(hic_sparse) * 1000 )
    ab_score=PC1$BPH456[bins$chr == chr][1:ncol(hic_sparse)]


    df <- data.frame(x=ab_score,y=lamin_contacts)
    r <- cor(df$x, df$y, use = "complete.obs",method = "spearman")

    gg <- ggplot(df, aes(x = x, y = y)) +
      geom_point() +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = sprintf("r = %.2f", r),  # two decimals
        hjust = 1.1, vjust = 1.5
      ) +
      theme_classic()+labs(x = "A/B score", y = "Lamin contacts")

    pdf(paste0("Hi-C/figures/AB-lamin-contact-scatter/BPH456-",chr,".pdf"), width = 6, height = 6)
    print(gg)
    dev.off()

}


