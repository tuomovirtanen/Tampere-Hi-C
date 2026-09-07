

library(anndata)
file="aggregated_sc_data_with_cell_type_annotation_20230613.h5ad"
adata <- read_h5ad(file)

cell_type = adata$obs[["refined_celltypes"]]
cluster=adata$obs[["VI_clusters"]]
sample=adata$obs[["sample"]]

crpc=unique(adata$obs[["sample"]][adata$obs[["phenotype"]] == "CRPC"])

table(cluster[cell_type=="Epithelial"])
keep_clusters=unique(cluster[cell_type=="Epithelial"])


basal_signature=read.table("signature-basal.txt")$V1
basal_signature=which(rownames(adata$var) %in% basal_signature)
luminal_signature=read.table("signature-luminal.txt")$V1
luminal_signature=which(rownames(adata$var) %in% luminal_signature)
cancer_signature=read.table("signature-cancer.txt")$V1
cancer_signature=which(rownames(adata$var) %in% cancer_signature)

norm_counts=adata$X
signature_scores <- data.frame(
  basal = Matrix::rowMeans(norm_counts[,basal_signature , drop = FALSE]),
  luminal = Matrix::rowMeans(norm_counts[,luminal_signature, drop = FALSE]),
  cancer = Matrix::rowMeans(norm_counts[,cancer_signature, drop = FALSE]),
  cluster = cluster,
  cell_type = cell_type
)


normals=adata$obs[["phenotype"]] == "normal"

normals[signature_scores$cell_type == "Epithelial"]
cluster[signature_scores$cell_type == "Epithelial"]
normals_by_cluster <- split(normals, cluster)
# calculate fraction of normal cells per cluster
percent_normal_per_cluster <- sapply(normals_by_cluster, function(x) sum(x) / length(x) * 100)
percent_normal_per_cluster

keep_clusters



write.table(percent_normal_per_cluster[names(percent_normal_per_cluster) %in% keep_clusters],"epithelial_clusters_normal_sample_fraction.txt")

#########
########

d <- read.table("epithelial_clusters_normal_sample_fraction.txt")

pdf("figures/single_cell_normal_cell_contribution.pdf")
barplot(
  d$x,
  names.arg = rownames(d),
  col = "steelblue",
  las = 2,                 # make x labels vertical
  ylab = "Fraction (%)",
  xlab = "Cluster",
  main = "Normal Sample Fraction by Epithelial Cluster"
)
dev.off()

signature_scores <- signature_scores[signature_scores$cell_type == "Epithelial", ]
library(ggplot2)
signature_scores <- read.table("epithelial_clusters_signature_scores.txt",header=TRUE)
signature_scores$cluster <- factor(signature_scores$cluster)

gg <- ggplot(signature_scores, aes(x = cluster, y = basal, fill = cluster)) +
  geom_violin(trim = FALSE, scale = "width") +
  theme_minimal() +
  theme(legend.position = "none")

pdf("figures/epithelial_clusters_signature_scores-basal.pdf", width = 10, height = 6)
print(gg)
dev.off()

gg <- ggplot(signature_scores, aes(x = cluster, y = luminal, fill = cluster)) +
  geom_violin(trim = FALSE, scale = "width") +
  theme_minimal() +
  theme(legend.position = "none")

pdf("figures/epithelial_clusters_signature_scores-luminal.pdf", width = 10, height = 6)
print(gg)
dev.off()


gg <- ggplot(signature_scores, aes(x = cluster, y = cancer, fill = cluster)) +
  geom_violin(trim = FALSE, scale = "width") +
  theme_minimal() +
  theme(legend.position = "none")

pdf("figures/epithelial_clusters_signature_scores-cancer.pdf", width = 10, height = 6)
print(gg)
dev.off()

#################

#
luminal_clusters <- c(5,12,13,27,28,37)
normal_luminal <- c(5)
cancer_luminal <- c(12,13,27,28,37)
basal_clusters <- c(9,10,16,23)

sum(cluster %in% normal_luminal & norm)
sum(cluster %in% cancer_luminal & !norm)
sum(cluster %in% basal_clusters & norm)

norm=sample %in% unique(sample[normals])

raw_counts=adata$layers[["counts"]]
out=list()
for(s in unique(sample[normals])){
  keep_s = sample == s & cluster %in% normal_luminal
  pseudo_bulk=Matrix::colSums(raw_counts[keep_s,])
  out[[s]] <- pseudo_bulk
}
out <- do.call(cbind, out)
write.table(out,"normal_luminal_pseudo_bulk.txt",sep="\t",row.names=TRUE,quote=FALSE)


out=list()
for(s in unique(sample[normals])){
  keep_s = sample == s & cluster %in% basal_clusters
  pseudo_bulk=Matrix::colSums(raw_counts[keep_s,])
  out[[s]] <- pseudo_bulk
}
out <- do.call(cbind, out)
write.table(out,"normal_basal_pseudo_bulk.txt",sep="\t",row.names=TRUE,quote=FALSE)



out=list()
for(s in unique(sample[!normals])){
  keep_s = sample == s & cluster %in% cancer_luminal
  if(sum(keep_s) < 2) next
  pseudo_bulk=Matrix::colSums(raw_counts[keep_s,])
  out[[s]] <- pseudo_bulk
}

out <- do.call(cbind, out)
write.table(out,"cancer_luminal_pseudo_bulk.txt",sep="\t",row.names=TRUE,quote=FALSE)






###############




library(edgeR)
normal_luminal <- read.table("normal_luminal_pseudo_bulk.txt",header=TRUE,row.names=1)
normal_luminal <- normal_luminal[,colSums(normal_luminal) > 0]
cancer_luminal <- read.table("cancer_luminal_pseudo_bulk.txt",header=TRUE,row.names=1)
cancer_luminal <- cancer_luminal[,colSums(cancer_luminal) > 0]
normal_basal <- read.table("normal_basal_pseudo_bulk.txt",header=TRUE,row.names=1)
normal_basal <- normal_basal[,colSums(normal_basal) > 0]

dataset_map=read.table("sample_dataset_map.txt",header=TRUE)
counts <- cbind(normal_luminal,normal_basal,cancer_luminal)

dataset <- c()
for(s in colnames(counts)){
  dataset <- c(dataset,dataset_map$dataset[match(s,dataset_map$sample)])
}

colnames(counts) <- c(paste0(colnames(normal_luminal),"_Normal"), paste0(colnames(normal_basal),"_Basal"), paste0(colnames(cancer_luminal),"_Cancer"))
group <- factor(c(rep("Normal", ncol(normal_luminal)),rep("Basal", ncol(normal_basal)), rep("Cancer", ncol(cancer_luminal))))
keep = dataset %in% c("chen_2022","hirz_2023","song_2022","wong_2022") & colSums(counts) > 50000
counts <- counts[,keep]
group <- group[keep]
dataset <- dataset[keep]

dge <- DGEList(counts = counts, group = group)
# Filter lowly expressed genes
#keep <- filterByExpr(dge, group = group, min.count = 1)
keep <- filterByExpr(dge, group = group, min.count = 1)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~0 + group + dataset)
colnames(design)[1:3] <- c(levels(group))
# Estimate dispersions
dge <- estimateDisp(dge,design)
# Fit the GLM
fit <- glmFit(dge,design)

contrast_normal_vs_cancer <- makeContrasts(Cancer - Normal, levels=design)
lrt_normal_vs_cancer <- glmLRT(fit, contrast=contrast_normal_vs_cancer)

contrast_basal_vs_normal <- makeContrasts(Normal - Basal, levels=design)
lrt_basal_vs_normal <- glmLRT(fit, contrast=contrast_basal_vs_normal)

# Extract top DE genes for each comparison
top_normal_vs_cancer <- topTags(lrt_normal_vs_cancer, n=Inf)$table
top_basal_vs_normal <- topTags(lrt_basal_vs_normal, n=Inf)$table

out1=data.frame(gene=rownames(top_normal_vs_cancer),logFC=top_normal_vs_cancer$logFC,FDR=top_normal_vs_cancer$FDR)
write.table(out1,"luminal_vs_cancer_pseudo_bulk_DE_results.txt",sep="\t",row.names=FALSE,quote=FALSE)

out2=data.frame(gene=rownames(top_basal_vs_normal),logFC=top_basal_vs_normal$logFC,FDR=top_basal_vs_normal$FDR)
write.table(out2,"basal_vs_luminal_pseudo_bulk_DE_results.txt",sep="\t",row.names=FALSE,quote=FALSE)