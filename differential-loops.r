library(edgeR)
library(limma)
library(pracma)
library(GenomicRanges)
library(locfit)
library(zoo)

path <- "HiC/loops/mustache-10kb"
inpath="HiC/loops/mustache-10kb/diff_loops"
outpath="HiC/loops/mustache-10kb/diff_loops/diff_loops_non_unique"


genome_size <- read.table("Hi-C/references/hg38/hg38.fa.chrom.sizes.sorted")
res=10000
loops <- read.table(paste0(path,"/loops/unified_loops.pgl"),header=T)
loop_quant <- read.table(paste0(inpath,"/loop_quant.txt"))

in_samples <- unlist(lapply(loops[,8],function(x) FUN = length(unlist(strsplit(x,",")))))

loops <- loops[in_samples > 1,]
loop_quant <- loop_quant[in_samples > 1,]

groups <- factor(c(rep("BPH",6),rep("CRPC",6),rep("PC",6)))

library(limma)
design = model.matrix(~0+groups)

colnames(design) = levels(groups)

myContrasts <- c("PC-BPH","CRPC-BPH","CRPC-PC")
contrast.matrix = eval(as.call(c(as.symbol("makeContrasts"),as.list(myContrasts),levels=list(design))))
contrast.matrix




coverage <- read.table(paste0(inpath,"/coverage.txt"))
coverage[,1] <- paste0("chr",coverage[,1])


cov_dge <- DGEList(counts=coverage[-c(1:3)],group=groups)
cov_dge <-  estimateDisp(cov_dge,design=design)

dge <- DGEList(counts=loop_quant,group=groups)
#dge$samples$lib.size <- lib.sizes[,1]
#dge$samples$lib.size <- cov_dge$samples$lib.size
#dge <- calcNormFactors(dge)
dge <-  estimateDisp(dge,design=design)


loop_up <-  makeGRangesFromDataFrame(loops,
                         keep.extra.columns=FALSE,
                         ignore.strand=TRUE,
                         seqnames.field="V1",
                         start.field="V2",
                         end.field="V3")
start(loop_up) <- start(loop_up) + res/2
end(loop_up) <- end(loop_up) - res/2

loop_down   <-  makeGRangesFromDataFrame(loops,
						 keep.extra.columns=FALSE,
                         ignore.strand=TRUE,
                         seqnames.field="V4",
                         start.field="V5",
                         end.field="V6")
start(loop_down) <- start(loop_down) + res/2
end(loop_down) <- end(loop_down) - res/2
						 
bins <-  makeGRangesFromDataFrame(coverage,
						 keep.extra.columns=FALSE,
                         ignore.strand=TRUE,
                         seqnames.field="chr",
                         start.field="start",
                         end.field="end")
						 
anchor1 <- as.matrix(findOverlaps(loop_up,bins))
anchor2	<- as.matrix(findOverlaps(loop_down,bins))

###this part is modified from diffHiC normalizeCNV
prior.count=3
span=0.3
maxk=500

cont.cor <- 0.5
cont.cor.scaled <- cont.cor * dge$samples$lib.size/mean(dge$samples$lib.size)

ab <- aveLogCPM(dge, prior.count=cont.cor)

mave <- aveLogCPM(cov_dge,prior.count=prior.count)
mave <- rollmean(mave,k=3)

mab <- cpm(cov_dge, log=TRUE, prior.count=prior.count)
#rolling mean because the loops were extended to 1 bin to each direction
mab <- apply(mab,2,function(x) FUN=rollmean(x,k=3))
mab <- mab - mave

ma.adjc <- mab[anchor1[,2],,drop=FALSE] 
mt.adjc <- mab[anchor2[,2],,drop=FALSE]

offsets <- matrix(0, nrow=nrow(dge), ncol=ncol(dge))

for(lib in 1:ncol(dge)){
	print(lib)
	ma.fc <- ma.adjc[,lib]
	mt.fc <- mt.adjc[,lib]
	mfc1 <- (ma.fc + mt.fc)/2
	mfc2 <- abs(ma.fc - mt.fc)
	all.cov <- list(mfc1, mfc2, ab)
	
	i.fc <- cpm(dge,log=T,prior.count=cont.cor)[,lib] - ab
	#i.fc <- log2(dge$counts[,lib] + cont.cor.scaled[lib]) - ab 
	cov.fun <- do.call(lp, c(all.cov, nn=span, deg=1))
	fit <- locfit(i.fc ~ cov.fun, maxk=maxk, lfproc=locfit.robust)
	
	offsets[,lib] <- fitted(fit)
	
	
	# ma.fc <- ma.adjc[,lib]
	# mt.fc <- mt.adjc[,lib]
	# mfc <- apply(cbind(ma.fc,mt.fc),1,function(x) FUN=x[which.min(abs(x))])

	# all.cov <- list(mfc, ab)
	# i.fc <- cpm(dge,log=T,prior.count=cont.cor)[,lib] - ab

	# cov.fun <- do.call(lp, c(all.cov, nn=span, deg=1))
	# fit <- locfit(i.fc ~ cov.fun, maxk=maxk, lfproc=locfit.robust)
	
	# png(paste0(outpath,"/plots/loc_fit",row.names(dge$samples)[lib],".png"),type="cairo")
	# par(mfrow=c(1,1))
	# plot(fit)
	# dev.off()
	
	# offsets[,lib] <- fitted(fit)
	
}

fits <- offsets

offsets <- 2^offsets
offsets <- log(offsets)
#offsets <- offsets - rowMeans(offsets)


cnv_dge <- DGEList(counts=loop_quant,group=groups)

o <- outer( rep(1,nrow(cnv_dge)), getOffset(cnv_dge)) + offsets


cnv_dge$offset <- o
cnv_dge <-  estimateDisp(cnv_dge,design=design)

cov <- apply(coverage[-c(1:3)],2,function(x) FUN=rollmean(x,k=3))
bins <- coverage[,1:3]
chrbreaks <- which(!!diff(as.numeric(as.factor(bins$chr)))) 
		
textbreaks <- c(0,chrbreaks)
textbreaks <- textbreaks[-length(textbreaks)] + diff(textbreaks) / 2


cp <- cpm(dge, log = TRUE, normalized.lib.sizes=TRUE)
ave_cp <- aveLogCPM(dge)

colnames(cp) <- colnames(loop_quant)
cp <- data.frame(cp)

cp <- cp-ave_cp


cp_cnv <- cpm(cnv_dge, log = TRUE)
ave_cp_cnv <- aveLogCPM(cnv_dge)

colnames(cp_cnv) <- colnames(loop_quant)
cp_cnv <- data.frame(cp_cnv)

cp_cnv <- cp_cnv-ave_cp_cnv

dir.create(file.path(outpath,"plots"), showWarnings = FALSE)
for(i in 1:ncol(dge)){

		s <- colnames(loop_quant)[i]

		logcov <- log2((cov[,i]+1)/median(cov[,i]))
		
		
		png(paste0(outpath,"/plots/",s,".png"),type="cairo",height=1200,width=3600)
		par(mfrow=c(4,1))
		plot(logcov,main=paste0(s," coverage"),ylim=c(-2,2.2),xlim=c(0,length(logcov)))
		abline(h=0)
		abline(v=chrbreaks)
		text(textbreaks,2.2, paste0("chr",c(1:22,"X","Y")),xpd=NA)

		cov_bias <- mab[,i]


		plot(anchor1[,2],cov_bias[anchor1[,2]],main=paste0(s," coverage"),ylim=c(-2,2.2),xlim=c(0,length(logcov)))
		abline(h=0)
		abline(v=chrbreaks)
		text(textbreaks,2.2, paste0("chr",c(1:22,"X","Y")),xpd=NA)
		
		
		plot(anchor1[,2],cp[,i],main=paste0(s, " log2f cpm"),ylab="norm cpm count",ylim=c(-2,2.2),xlim=c(0,length(logcov)))
		abline(h=0)
		abline(v=chrbreaks)
		text(textbreaks,2.2, paste0("chr",c(1:22,"X","Y")),xpd=NA)
		points(anchor1[,2],fits[,i],col="red")

		plot(anchor1[,2],cp_cnv[,i],main=paste0(s, " log2f cpm"),ylab="norm cpm count",ylim=c(-2,2.2),xlim=c(0,length(logcov)))
		abline(h=0)
		abline(v=chrbreaks)
		text(textbreaks,2.2, paste0("chr",c(1:22,"X","Y")),xpd=NA)
		dev.off()

		# plot(anchor1[,2],cp[,i]-fits[,i],main=paste0(s, " log2f cpm"),ylab="norm cpm count",ylim=c(-2,2.2),xlim=c(0,length(logcov)))
		# abline(h=0)
		# abline(v=chrbreaks)
		# text(textbreaks,2.2, paste0("chr",c(1:22,"X","Y")),xpd=NA)
		# dev.off()

}

cp <- cpm(cnv_dge, log = TRUE, normalized.lib.sizes=TRUE)
colnames(cp) <- colnames(loop_quant)
write.table(cp,paste0(outpath,"/edgeR_cpm_norm.txt"))

dgefit <- glmQLFit(cnv_dge,design=design)
results <- list()
dge_results <- list()
for(i in 1:ncol(contrast.matrix)){

	dge_results[[i]] <- glmQLFTest(dgefit,contrast=contrast.matrix[,i])
	results[[i]] <- dge_results[[i]]$table
	results[[i]]$fdr <- p.adjust(results[[i]]$PValue,method="fdr")
	colnames(results[[i]]) <- paste0(colnames(contrast.matrix)[i],"-",colnames(results[[i]]))
}

results <- do.call(cbind,results)

write.table(results,paste0(outpath,"/edgeR_results.txt"))
write.table(offsets,paste0(outpath,"/edgeR_cnv_offsets.txt"))
