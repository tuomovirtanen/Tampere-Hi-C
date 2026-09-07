
AB_scores <- read.table("PC1.txt")
outdir <- "limma"

library(limma)
library(gplots)

b <- AB_scores[,-c(1:3)]*100
keep=rowSums(is.na(b)) == 0
b<-b[keep,]
row.names(b) <- which(keep)
groups <- factor(c(rep("BPH",6),rep("CRPC",6),rep("PC",6)))

BPH <- groups == "BPH"
PC <- groups == "PC"
CRPC <- groups == "CRPC"

cf <- c(
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	0.798,
	0.645,
	0.575,
	0.279,
	0.366,
	0.469,
	0.231,
	0.322,
	0.437,
	0.706,
	0.163,# estimate based on  TMPRSS-ERG fusion deletion segment
	0.390 # estimate based on  TMPRSS-ERG fusion deletion segment
)

library(limma)
design = model.matrix(~0+groups)

colnames(design) = c(levels(groups))
design[7:18,1] = design[7:18,1] + (1-cf[7:18]) 
design[,2] = design[,2]*cf
design[,3] = design[,3]*cf

myContrasts <- c("PC-BPH","CRPC-BPH","CRPC-PC")
contrast.matrix = eval(as.call(c(as.symbol("makeContrasts"),as.list(myContrasts),levels=list(design))))
contrast.matrix=cbind(contrast.matrix,c(-2,1,1))
colnames(contrast.matrix)[4] <- "Cancer-BPH"
contrast.matrix

fit = lmFit(b, design)

fit2 = contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2,trend = F)
results <- decideTests(fit2)

results_p <- decideTests(fit2,adjust.method = "none",p.value = 0.01)

res1 <- topTable(fit2,coef="PC-BPH",number=nrow(b))
res1 <- res1[order(as.numeric(row.names(res1))),]
res2 <- topTable(fit2,coef="CRPC-BPH",number=nrow(b))
res2 <- res2[order(as.numeric(row.names(res2))),]
res3 <- topTable(fit2,coef="CRPC-PC",number=nrow(b))
res3 <- res3[order(as.numeric(row.names(res3))),]
res4 <- topTable(fit2,coef="Cancer-BPH",number=nrow(b))
res4 <- res4[order(as.numeric(row.names(res4))),]


out <- cbind(res1$P.Value,res2$P.Value,res3$P.Value,res4$P.Value,
			res1$adj.P.Val,res2$adj.P.Val,res3$adj.P.Val,res4$adj.P.Val,
			res1$logFC,res2$logFC,res3$logFC,res4$logFC)
colnames(out) <- c("p.val-PC-BPH","p.val-CRPC-BPH","p.val-CRPC-PC","p.val-Cancer-BPH",
"q.val-PC-BPH","q.val-CRPC-BPH","q.val-CRPC-PC","q.val-Cancer-BPH",
"log2f-PC-BPH","log2f-CRPC-BPH","log2f-CRPC-PC","log2f-Cancer-BPH")

write.table(out,paste0(outdir,"/limma_results.txt"))
write.table(res1,paste0(outdir,"/PC_BPH_limma_results.txt"))
write.table(res2,paste0(outdir,"/CRPC_BPH_limma_results.txt"))
write.table(res3,paste0(outdir,"/CRPC_PC_limma_results.txt"))
write.table(res4,paste0(outdir,"/Cancer-BPH_limma_results.txt"))
