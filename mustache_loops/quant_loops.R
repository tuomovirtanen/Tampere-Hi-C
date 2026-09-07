
library(Matrix)

path <- "HiC/loops/mustache-10kb"
outpath="HiC/loops/mustache-10kb/diff_loops"
hic_files <- "/hic_files/v_14/"
res=10000


loops <- read.table(paste0(path,"/loops/unified_loops.pgl"),stringsAsFactors=F,header=T)


samples <- unique(unlist(strsplit(loops[,8],",")))
samples <- samples[order(samples)]
chrs <- c(1:22,"X","Y")

loop_count <- matrix(NA,ncol=length(samples),nrow=nrow(loops))
coverage <- list()
lib_size <- rep(list(0),length(samples))
names(lib_size) <- samples

for(chr in chrs){
	print(chr)
	chr_coverage <- list()
	chr_loops <- loops[,1] == paste0("chr",chr)
	
		loop_idx <- (loops[chr_loops,c(2,3,5,6)]/res)+1
		loop_idx[,c(2,4)] <- loop_idx[,c(2,4)]-1
		
		anch1_pad <- loop_idx[,2]-loop_idx[1] < 2
		
		anch2_pad <- loop_idx[,4]-loop_idx[3] < 2 
		
		loop_idx[anch1_pad,1] <- loop_idx[anch1_pad,1] - 1
		loop_idx[anch1_pad,2] <- loop_idx[anch1_pad,2] + 1
		loop_idx[anch2_pad,3] <- loop_idx[anch2_pad,3] - 1
		loop_idx[anch2_pad,4] <- loop_idx[anch2_pad,4] + 1
		
	i=0	
	for(sample in samples){
		print(sample)
		i=i+1
		mat <- strawr::straw("observed","NONE",paste0(hic_files,"/",sample,".hic"),chr,chr,"BP",res)
		lib_size[[sample]] <- lib_size[[sample]] + sum(mat[,3],na.rm=TRUE)
		mat[,1:2] <- mat[,1:2]/res+1
		d <- max(max(mat[,1]),max(mat[,2]))
		mat <- sparseMatrix(i=mat[,1],j=mat[,2],x=mat[,3],dims=c(d,d),symmetric=TRUE)
		mat <- as.matrix(mat)

		loop_count[chr_loops,i] <- apply(loop_idx,1,function(x) FUN=sum(mat[(x[[1]]:x[[2]]),(x[[3]]:x[[4]])]))
		
		chr_coverage[[sample]] <- rowSums(mat)
	}
	
	lim <- min(unlist(lapply(chr_coverage,length)))
	chr_coverage <- lapply(chr_coverage,function(x) FUN=x[1:lim])
	chr_coverage <- do.call(cbind,chr_coverage)
	
	bins <- data.frame(chr=chr,start=seq(0,by=res,length.out=lim),end=seq(res,by=res,length.out=lim))
	coverage[[chr]] <- cbind(bins,chr_coverage)
}

coverage <- do.call(rbind,coverage)
colnames(loop_count) <- samples

write.table(loop_count,paste0(outpath,"/loop_quant.txt"))
write.table(coverage,paste0(outpath,"/coverage.txt"))

libsize <- unlist(lib_size)
write.table(libsize,paste0(outpath,"/library_size.txt"))

