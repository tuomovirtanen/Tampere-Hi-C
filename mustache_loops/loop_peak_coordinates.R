args = commandArgs(trailingOnly = TRUE)

unified<-args[1]
sorted <- args[2]
out <- args[3]
res <- as.numeric(args[4])

options(scipen=100)
unified <- read.table(unified)
sorted <- read.table(sorted)




starts <- rep(NA,nrow(unified))
ends <- rep(NA,nrow(unified))
for(i in 1:nrow(unified)){
	idx <- as.numeric(unlist(strsplit(as.character(unified[i,10]),",")))
	starts[i] <- round(median(sorted[idx,2])/res)*res
	ends[i] <- round(median(sorted[idx,5])/res)*res
}


unified[,2] <- starts
unified[,3] <- starts + res
unified[,5] <- ends
unified[,6] <- ends + res
write.table(unified,out,sep="\t",quote=F,row.names=FALSE,col.names=TRUE)