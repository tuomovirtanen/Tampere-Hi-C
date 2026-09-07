library(doParallel)
library(Matrix)

args = commandArgs(trailingOnly = TRUE)
#(input,chr,res,matrixOut)
hic_file <-args[1]
segments_file <- args[2]
output <- args[3]

res=100000

chrs <- c(1:22,"X","Y")
#chrs <- as.character(c(20:22))
data=list()
for(chr in chrs){
	mat <- strawr::straw("SCALE",hic_file,as.character(chr),as.character(chr),"BP",res,"observed")
	mat[,1:2] <- mat[,1:2]/res+1
	mat_oe <- strawr::straw("SCALE",hic_file,as.character(chr),as.character(chr),"BP",res,"oe")
	data[[chr]] <- list(mat,mat_oe)
}

names(data)

out <- lapply(names(data), function(x) cbind(data[[x]][[1]],x))
names(out) <- names(data)
mat <- lapply(data,function(x)x[[1]])
oe <- lapply(data,function(x)x[[2]])

# if(!all(dim(mat) == dim(oe))){
# 	stop("observed and Oe matrix size differ")
# }

segments <- read.table(segments_file,header=T)

segments <- segments[,1:3]
segments[,1] <- gsub("chr","",segments[,1])
#For segment calls bin was set to center hence - res/2
segments[,2:3] <- ((segments[,2:3] - res/2) / res) +1


seg1=sapply(1:nrow(segments),function(x) paste0(segments[x,1],"-",segments[x,2]))
seg2=sapply(1:nrow(segments),function(x) paste0(segments[x,1],"-",segments[x,3]))

seg=unique(c(seg1,seg2))

seg_chrs=unlist(lapply(seg,function(x) strsplit(x,"-")[[1]][1]))
seg_bins=unlist(lapply(seg,function(x) as.numeric(strsplit(x,"-")[[1]][2])))

ord=order(seg_bins)
seg_chrs=seg_chrs[ord]
seg_bins=seg_bins[ord]
ord=order(seg_chrs)
seg_chrs=seg_chrs[ord]
seg_bins=seg_bins[ord]



segment_list <- list()
norms <- list()
k=1
for(i in 1:(length(seg)-1)){
	for(j in 1:(length(seg)-1)){
		#print(k)
		s11=seg_bins[i]
		s12=seg_bins[i+1]
		s21=seg_bins[j]
		s22=seg_bins[j+1]
		
		chr1=seg_chrs[i]
		chr2=seg_chrs[i+1]
		chr3=seg_chrs[j]
		chr4=seg_chrs[j+1]
	
		if(all(c(chr1,chr2,chr3)==chr4) & abs(s12-s11) > 2 & abs(s22-s21) > 2){
				pick = 	mat[[chr1]][,1] >= s11 & mat[[chr1]][,1] <= s12 & mat[[chr1]][,2] >= s21 & mat[[chr1]][,2] <= s22
				norm <- mean(oe[[chr1]][pick,3],na.rm=T)
				out[[chr1]][pick,3] <- mat[[chr1]][pick,3]/norm
				norms[[k]] <- norm
				k=k+1

		}
	}	
}

mat_normed=do.call(rbind,out)
juicer_out_block_normed <- data.frame(str1=0,
							chr1=mat_normed[,4],
							pos1=(mat_normed[,1]-1)*res+res/2,
							frag1=0,
							str2=1,
							chr2=mat_normed[,4],
							pos2=(mat_normed[,2]-1)*res+res/2,
							frag2=1,
							score=mat_normed[,3])

#juicer_out_block_normed
write.table(juicer_out_block_normed,output,quote=F,row.names=F,col.names=F,sep="\t")



