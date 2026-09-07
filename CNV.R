

library(doParallel)
library(Matrix)
library(edgeR)
library(copynumber)
hic_files <- "HiC/hic_files/"
res=100000
output <- "HiC/CNV/100kb"


samples <- dir(hic_files)
samples <- gsub(".hic","",samples)
samples <- samples[order(samples)]
chrs <- c(1:22,"X","Y")
#chrs <- as.character(c(20:22))
coverage <- list()

no_cores <- min(detectCores() - 1,8) 

cl <- makeCluster(no_cores)  
registerDoParallel(cl)

cov_mat <- foreach(chr=chrs, .packages=c("Matrix")) %dopar% {
	chr_coverage <- list()

	i=0	
	for(sample in samples){
		print(sample)
		i=i+1
		mat <- strawr::straw("observed","NONE",paste0(hic_files,"/",sample,".hic"),as.character(chr),as.character(chr),"BP",res)
		mat[,1:2] <- mat[,1:2]/res+1
		d <- max(max(mat[,1]),max(mat[,2]))
		mat <- sparseMatrix(i=mat[,1],j=mat[,2],x=mat[,3],dims=c(d,d),symmetric=TRUE)
		mat <- as.matrix(mat)
		chr_coverage[[sample]] <- rowSums(mat)
		
	}

	lim <- min(unlist(lapply(chr_coverage,length)))
	chr_coverage <- lapply(chr_coverage,function(x) FUN=x[1:lim])
	chr_coverage <- do.call(cbind,chr_coverage)
	
	bins <- data.frame(chr=chr,start=seq(0,by=res,length.out=lim),end=seq(res,by=res,length.out=lim))
	coverage <- cbind(bins,chr_coverage)
	
	
	coverage
}
cov_mat <- do.call(rbind,cov_mat)

mat <- cov_mat[,-c(1:3)]
filt <- t(t(mat[,1:6])/colMeans(mat[,1:6]))
filt <- rowSums(filt < 0.1) == 0


mat <- mat[filt,]

mat_norm <- t(t(mat) / colSums(mat))
#mat_norm <- mat_norm / rowMeans(mat_norm[,1:6])
mat_norm <- mat_norm / apply(mat_norm[,1:6],1,median)
mat_norm <- log2(mat_norm)

mat_norm[mat_norm < -2] <- -2


idx <- as.numeric(row.names(mat_norm))
i=1
for(s in colnames(mat_norm)){
	png(paste0(s,".png"),type="cairo",width=2400,height=400)
	plot(idx,mat_norm[,i],main=s)
	dev.off()
	i=i+1
}

################
gamma=1000

i=0
for(sample in samples){
	i=i+1
	d <- data.frame(Chrom=cov_mat[filt,1],pos=cov_mat[filt,2],x=mat_norm[,i])

	d <- d[!is.na(d$x),]
	d <- d[d$x > -2,]

	d.wins <- winsorize(data=d,verbose=FALSE)

	seg.wins <- pcf(data=d.wins,gamma=gamma)

	png(paste0(output,"/",sample,"_segments.png"),height = 1200,width=1200,type="cairo")
	plotSample(data=d.wins,segments=seg.wins,layout=c(5,5),sample=1,cex=3)
	dev.off()

	cn <- rep(NA,nrow(cov_mat))

	for(j in 1:nrow(seg.wins)){
	  start = which(cov_mat[,1]==seg.wins[j,2] & cov_mat[,2]==seg.wins[j,4])
	  end = which(cov_mat[,1]==seg.wins[j,2] & cov_mat[,2]==seg.wins[j,5])
	  cn[start:end] <- seg.wins[j,7]
	}


	write.table(seg.wins,paste0(output,"/",sample,"_segments.txt"))
	cn_log2=cn
	cn=(2^cn)*2
	cn[is.na(cn)] <- 0
	cn_log2[is.na(cn_log2)] <- 0
	out <- cbind(cov_mat[,1:3],cn,cn_log2,ifelse(filt,1,0))
	colnames(out) <- c("chr","start","end","cnv","cnv_log2","filt")
	write.table(out,paste0(output,"/",sample,"_cn.txt"),row.names = F,
				col.names = T,quote = F,sep="\t")
}





