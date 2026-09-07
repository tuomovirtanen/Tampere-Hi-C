library(GenomicRanges)


SVs=read.table("Hi-C/data/2018_04_15_list_manta_SV.txt",header = T,stringsAsFactors = F)
compartments=read.table("Hi-C/Compartment_classes_Cancer-BPH.bed",header=TRUE)
coords=compartments[,1:3]


AB = rep(NA,nrow(compartments))
AB[grepl("_A",compartments$class)] = "A"
AB[grepl("_B",compartments$class)] = "B"


b_gr = makeGRangesFromDataFrame(coords,
                                 keep.extra.columns=FALSE,
                                 ignore.strand=TRUE,
                                 seqnames.field="chr",
                                 start.field="start",
                                 end.field="end")
seqlevelsStyle(b_gr) <- "UCSC"


sv=SVs[SVs$svtype == "BND",]

sv_breaks=list()
for(i in 1:nrow(sv)){
  if( sv$chrom_start[i] < sv$chrom_end[i]){
    sv_breaks[[i]] <- unlist(sv[i,5:8])
  } else {
    sv_breaks[[i]] <- unlist(sv[i,c(7,8,5,6)])
  }
}


sv_breaks = do.call(rbind,sv_breaks)
sv_breaks=unique(sv_breaks)
sv_breaks=data.frame(sv_breaks)
sv_gr1 = makeGRangesFromDataFrame(sv_breaks,
                                 keep.extra.columns=FALSE,
                                 ignore.strand=TRUE,
                                 seqnames.field="chrom_start",
                                 start.field="pos_start",
                                 end.field="pos_start")


sv_gr2 = makeGRangesFromDataFrame(sv_breaks,
                                  keep.extra.columns=FALSE,
                                  ignore.strand=TRUE,
                                  seqnames.field="chrom_end",
                                  start.field="pos_end",
                                  end.field="pos_end")



sv_ab1 <- rep(NA,length(sv_gr1))
ol1 = findOverlaps(sv_gr1,b_gr)
sv_ab1[queryHits(ol1)] = AB[subjectHits(ol1)]

sv_ab2 <- rep(NA,length(sv_gr2))
ol2 = findOverlaps(sv_gr2,b_gr)
sv_ab2[queryHits(ol2)] = AB[subjectHits(ol2)]

tb=table(sv_ab1,sv_ab2)
tb


n=1000
exp_a=0
exp_b=0
exp_ab=0
a_p=0
b_p=0
ab_p=0

rndm1_mat <- list()
rndm2_mat <- list()

for(i in which(!is.na(sv_ab1) & !is.na(sv_ab2))){

  chr1=sv_breaks$chrom_start[i]
  chr2=sv_breaks$chrom_end[i]
  
  rndm1=sample(AB[!is.na(AB) & coords$chr == chr1],n,replace=T)
  rndm2=sample(AB[!is.na(AB) & coords$chr == chr2],n,replace=T)

  rndm1_mat[[i]] <- rndm1
  rndm2_mat[[i]] <- rndm2
}

rndm1_mat <- do.call(rbind,rndm1_mat)
rndm2_mat <- do.call(rbind,rndm2_mat)


for(i in 1:ncol(rndm1_mat)){
  rndm1=rndm1_mat[,i]
  rndm2=rndm2_mat[,i]
  rtb=table(rndm1,rndm2)

  exp_a = exp_a + rtb[1,1]
  exp_b = exp_b + rtb[2,2]
  exp_ab = exp_ab + rtb[1,2] +rtb[2,1]

  a_p = a_p + ifelse((rtb[1,1] / sum(rtb)) > (tb[1,1] / sum(tb)),1,0)
  b_p = b_p + ifelse((rtb[2,2] / sum(rtb)) > (tb[2,2] / sum(tb)),1,0)
  ab_p = ab_p + ifelse(((rtb[1,2] + rtb[2,1]) / sum(rtb)) > ((tb[1,2] + tb[2,1]) / sum(tb)),1,0)
}


exp_a=exp_a / n
exp_b=exp_b / n
exp_ab=exp_ab /n
a_p=a_p / n
b_p=b_p / n
ab_p=ab_p / n


library(ggplot2)
gg_data=data.frame(class=factor(c("AA","BB","AB"),levels=c("AA","BB","AB")),count=c(tb[1,1]/exp_a,tb[2,2]/exp_b,(tb[1,2]+tb[2,1])/exp_ab))
gg <- ggplot()+geom_col(data=gg_data,aes(x=class,y=count))+
theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
ylab("observed / expected")

pdf("Hi-C/figures/SVs-AB.pdf")
print(gg)
dev.off()


#############

library(GenomicRanges)
library(Matrix)

SVs=read.table("Hi-C/data/2018_04_15_list_manta_SV.txt",header = T,stringsAsFactors = F)
compartments=read.table("Hi-C/Compartment_classes_Cancer-BPH.bed",header=TRUE)
coords=compartments[,1:3]

chroms <- c(paste0("chr", 1:22))
chroms

AB = rep(NA,nrow(compartments))
AB[grepl("_A",compartments$class)] = "A"
AB[grepl("_B",compartments$class)] = "B"


b_gr = makeGRangesFromDataFrame(coords,
                                 keep.extra.columns=FALSE,
                                 ignore.strand=TRUE,
                                 seqnames.field="chr",
                                 start.field="start",
                                 end.field="end")
seqlevelsStyle(b_gr) <- "UCSC"


sv=SVs[SVs$svtype == "BND",]
sv_breaks=list()
for(i in 1:nrow(sv)){
  if( sv$chrom_start[i] < sv$chrom_end[i]){
    sv_breaks[[i]] <- unlist(sv[i,5:8])
  } else {
    sv_breaks[[i]] <- unlist(sv[i,c(7,8,5,6)])
  }
}

head(sv_breaks)
sv_breaks = do.call(rbind,sv_breaks)
sv_breaks=data.frame(unique(sv_breaks))
sv_breaks=sv_breaks[sv_breaks$chrom_start %in% chroms & sv_breaks$chrom_end %in% chroms,]
sv_breaks$pos_start <- as.numeric(sv_breaks$pos_start)
sv_breaks$pos_end <- as.numeric(sv_breaks$pos_end)


sv_gr1 = makeGRangesFromDataFrame(sv_breaks,
                                 keep.extra.columns=FALSE,
                                 ignore.strand=TRUE,
                                 seqnames.field="chrom_start",
                                 start.field="pos_start",
                                 end.field="pos_start")


sv_gr2 = makeGRangesFromDataFrame(sv_breaks,
                                  keep.extra.columns=FALSE,
                                  ignore.strand=TRUE,
                                  seqnames.field="chrom_end",
                                  start.field="pos_end",
                                  end.field="pos_end")

sv_ab1 <- rep(NA,length(sv_gr1))
ol1 = findOverlaps(sv_gr1,b_gr)
sv_ab1[queryHits(ol1)] = AB[subjectHits(ol1)]

sv_ab2 <- rep(NA,length(sv_gr2))
ol2 = findOverlaps(sv_gr2,b_gr)
sv_ab2[queryHits(ol2)] = AB[subjectHits(ol2)]


tb=table(sv_ab1,sv_ab2)


chroms

tmp1=sv_breaks[, 1:2]
colnames(tmp1) <- c("chr","pos")
tmp2=sv_breaks[, 3:4]
colnames(tmp2) <- c("chr","pos")

breaks <- rbind(
  tmp1,
  tmp2
)

rndm_breaks_1 <- list()
rndm_breaks_2 <- list()
rnd_pair=rep(NA, nrow(sv_breaks))
for(i in 1:nrow(sv_breaks)){
  brk1=sample(1:nrow(breaks),1)
  chr1=breaks$chr[brk1]
  brk2=sample(which(breaks$chr != chr1),1)

  rndm_breaks_1[[i]] <- breaks[brk1,]
  rndm_breaks_2[[i]] <- breaks[brk2,]
}
rndm_breaks_1 <- do.call(rbind,rndm_breaks_1)
rndm_breaks_2 <- do.call(rbind,rndm_breaks_2)

#sv_breaks_random <- cbind(sv_breaks[,1:2],sv_breaks[rnd_pair,3:4])
sv_breaks_random <- cbind(rndm_breaks_1,rndm_breaks_2)
colnames(sv_breaks_random) <- colnames(sv_breaks)

rndm_sv_gr1 = makeGRangesFromDataFrame(sv_breaks_random,
                                 keep.extra.columns=FALSE,
                                 ignore.strand=TRUE,
                                 seqnames.field="chrom_start",
                                 start.field="pos_start",
                                 end.field="pos_start")


rndm_sv_gr2 = makeGRangesFromDataFrame(sv_breaks_random,
                                  keep.extra.columns=FALSE,
                                  ignore.strand=TRUE,
                                  seqnames.field="chrom_end",
                                  start.field="pos_end",
                                  end.field="pos_end")

rndm_sv_ab1 <- rep(NA,length(rndm_sv_gr1))
ol1 = findOverlaps(rndm_sv_gr1,b_gr)
rndm_sv_ab1[queryHits(ol1)] = AB[subjectHits(ol1)]

rndm_sv_ab2 <- rep(NA,length(rndm_sv_gr2))
ol2 = findOverlaps(rndm_sv_gr2,b_gr)
rndm_sv_ab2[queryHits(ol2)] = AB[subjectHits(ol2)]





hic_files <- list.files("Hi-C/hic_files/",pattern = "*.hic",full.names = T)

all_contacts=list()
all_contacts_rndm=list()

for(hic_file in hic_files){
  print(hic_file)
  contacts=rep(NA,length(sv_gr1))
  contacts_rndm=rep(NA,length(rndm_sv_gr1))
  sample_name = gsub(".hic","",basename(hic_file))
  for(i in 1:(length(chroms)-1)){
    for(j in (i+1):length(chroms)){
      chr1=chroms[i]
      chr2=chroms[j]

      print(chr1)
      print(chr2)

      mat <- strawr::straw("SCALE",hic_file,gsub("chr","",chr1),gsub("chr","",chr2),"BP",1000000,"oe")
      mat$x <- floor(mat$x / 1000000) +1
      mat$y <- floor(mat$y / 1000000) +1
      hic_sparse <- sparseMatrix(
        i = mat$x,
        j = mat$y,
        x = mat$counts,
        dims = c(max(mat$x), max(mat$y)),
        symmetric = FALSE # or TRUE for cis-chromosomal maps if appropriate
      )

      idx=which((sv_breaks$chrom_start == chr1 & sv_breaks$chrom_end == chr2) | (sv_breaks$chrom_start == chr2 & sv_breaks$chrom_end == chr1))
      for(k in idx){
        x_idx=floor(ifelse(sv_breaks$chrom_start[k] == chr1, sv_breaks$pos_start[k], sv_breaks$pos_end[k]) / 1000000) +1
        y_idx=floor(ifelse(sv_breaks$chrom_start[k] == chr1, sv_breaks$pos_end[k],sv_breaks$pos_start[k]) / 1000000) +1
        contacts[k] <- hic_sparse[x_idx,y_idx]
      }

      idx=which((sv_breaks_random$chrom_start == chr1 & sv_breaks_random$chrom_end == chr2) | (sv_breaks_random$chrom_start == chr2 & sv_breaks_random$chrom_end == chr1))
      for(k in idx){
        x_idx=floor(ifelse(sv_breaks_random$chrom_start[k] == chr1, sv_breaks_random$pos_start[k], sv_breaks_random$pos_end[k]) / 1000000) +1
        y_idx=floor(ifelse(sv_breaks_random$chrom_start[k] == chr1, sv_breaks_random$pos_end[k],sv_breaks_random$pos_start[k]) / 1000000) +1
        contacts_rndm[k] <- hic_sparse[x_idx,y_idx]
      }
    }
  }
  class=rep(NA,length(sv_gr1))
  class[sv_ab1 == "A" & sv_ab2 == "A"] <- "AA"
  class[sv_ab1 == "B" & sv_ab2 == "B"] <- "BB"
  class[(sv_ab1 == "A" & sv_ab2 == "B") | (sv_ab1 == "B" & sv_ab2 == "A")] <- "AB"


  class_rndm=rep(NA,length(sv_gr1))
  class_rndm[rndm_sv_ab1 == "A" & rndm_sv_ab2 == "A"] <- "rand-AA"
  class_rndm[rndm_sv_ab1 == "B" & rndm_sv_ab2 == "B"] <- "rand-BB"
  class_rndm[(rndm_sv_ab1 == "A" & rndm_sv_ab2 == "B") | (rndm_sv_ab1 == "B" & rndm_sv_ab2 == "A")] <- "rand-AB"

  head(contacts)
  head(class)

  library(ggplot2)
  df <- data.frame(
    class = factor(c(class,class_rndm),levels=c("AA" ,"rand-AA","BB" ,"rand-BB","AB" ,"rand-AB")),   # make sure class is a factor
    contacts = c(contacts,contacts_rndm)
  )
  gg <- ggplot(df, aes(x = class, y = contacts)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw()+coord_cartesian(ylim=c(0,4))+

  pdf(paste0("Hi-C/figures/compartment-sv-distance-",sample_name,".pdf"))
  print(gg)
  dev.off()

  all_contacts[[sample_name]] <- contacts
  all_contacts_rndm[[sample_name]] <- contacts_rndm
}

avg_contacts = rowMeans(do.call(cbind, all_contacts),na.rm=T)
avg_contacts_rndm = rowMeans(do.call(cbind, all_contacts_rndm),na.rm=T)

df <- data.frame(
  class = factor(c(class,class_rndm),levels=c("AA" ,"rand-AA","BB" ,"rand-BB","AB" ,"rand-AB")),   # make sure class is a factor
  contacts = c(avg_contacts,avg_contacts_rndm)
)
gg <- ggplot(df, aes(x = class, y = contacts)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw()+coord_cartesian(ylim=c(0,4))

pdf(paste0("Hi-C/figures/compartment-sv-distance-average.pdf"))
print(gg)
dev.off()