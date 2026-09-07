
library(GenomicRanges)
compartments=read.table("Compartment_classes_Cancer-BPH.bed",header=TRUE)

shuffle_blocks <- function(x) {
  # Perform run-length encoding to identify blocks of consecutive, identical values.
  # rle() returns a list with two components:
  #   - lengths: a vector of the lengths of the runs
  #   - values: a vector of the corresponding values
  rle_result <- rle(x)
  # Get the number of blocks to be shuffled.
  num_blocks <- length(rle_result$values)
  # Create a vector of shuffled indices from 1 to the number of blocks.
  shuffled_indices <- sample(1:num_blocks, num_blocks)
  # Reorder the values and lengths vectors using the shuffled indices.
  shuffled_rle <- list(
    values = rle_result$values[shuffled_indices],
    lengths = rle_result$lengths[shuffled_indices]
  )
  # Use inverse run-length encoding to reconstruct the vector from the shuffled
  # blocks and their corresponding lengths.
  result_vector <- inverse.rle(shuffled_rle)
  # Return the final shuffled vector.
  return(result_vector)
}

atac_tf <- read.table("new-TF-anno-consensus_peaks_cleaned_CPM_cutoff_4.51_n_85450.bed",header=TRUE)

hoxb13_columns = grepl("\\.HOXB13\\.",colnames(atac_tf))
ar_columns = grepl("\\.AR\\.",colnames(atac_tf))
foxa1_columns = grepl("\\.FOXA1\\.",colnames(atac_tf))

hoxb13=rowSums(atac_tf[,hoxb13_columns]) > 0.25*sum(hoxb13_columns)
ar=rowSums(atac_tf[,ar_columns]) > 0.25*sum(ar_columns)
foxa1=rowSums(atac_tf[,foxa1_columns]) > 0.25*sum(foxa1_columns)

comp_list <- split(compartments[,4], compartments[,1])
ATAC_peaks <- read.table("consensus_peaks_cleaned_CPM_cutoff_4.51_n_85450.bed",header=F)
colnames(ATAC_peaks) =gsub("_","",colnames(ATAC_peaks))
peaks_bed <- ATAC_peaks[,1:3]
atac_gr <- makeGRangesFromDataFrame(peaks_bed,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=TRUE,
                                    seqnames.field="V1",
                                    start.field="V2",
                                    end.field="V3")
seqlevelsStyle(atac_gr)<- "UCSC"

bins_gr <- makeGRangesFromDataFrame(compartments,
                                   keep.extra.columns=FALSE,
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field="chr",
                                   start.field="start",
                                   end.field="end")
seqlevelsStyle(bins_gr)<- "UCSC"



ols=atac_gr %over% bins_gr[compartments[,4] == "Activating_B"]

#############
anno_files <- list.files("/data/chromHMM",full.names=TRUE)
samples <- basename(anno_files);samples <- gsub("_ChromHMM_dense_regressed_sorted_1.bed.processed.bed","",samples)
anno_list <- list()
for(k in 1:4){
  annotations <- read.table(anno_files[k])
  anno_gr <- makeGRangesFromDataFrame(annotations,
                                        keep.extra.columns=FALSE,
                                        ignore.strand=TRUE,
                                        seqnames.field="V1",
                                        start.field="V2",
                                        end.field="V3")
  anno <- rep(NA,length(atac_gr))
  ol <- findOverlaps(atac_gr,anno_gr)
  anno[queryHits(ol)] <- annotations[subjectHits(ol),4]

  anno_list[[k]] <- anno
}

anno_list <- do.call(cbind,anno_list)
anno_list[is.na(anno_list)] <- "None"
cancer_enh <- (anno_list[,1] == "Enhancer" | anno_list[,4] == "Enhancer") & (anno_list[,3] != "Enhancer") 
Normal_enh <- (anno_list[,3] == "Enhancer") & !(anno_list[,1] == "Enhancer" | anno_list[,4] == "Enhancer")

##############


data <- list(hoxb13,ar,foxa1,cancer_enh,Normal_enh)
results <- data.frame(
  x = character(),
  odds = numeric(),
  p = numeric(),
  stringsAsFactors = FALSE
)

k=0
for(x in data){
  k <- k +1
  print(k)
  tb = table(x, ols)
  odds = log2((tb[1]/tb[2]) / (tb[3]/tb[4]))

  rndm_odds = rep(NA, 10000)
  for(i in 1:10000){
    shuffeled = unlist(lapply(comp_list, function(x) shuffle_blocks(x)))
    tb_rndm = table(x, atac_gr %over% bins_gr[shuffeled == "Activating_B"])
    rndm_odds[i] = log2((tb_rndm[1]/tb_rndm[2]) / (tb_rndm[3]/tb_rndm[4]))
  }

  p = sum(abs(rndm_odds) >= abs(odds))

  results <- rbind(results,
                   data.frame(x = k,
                              odds = odds,
                              p = p))
  print(odds)
  print(p)
}

write.table(results, "compartment_permutation_results.txt",
            sep = "\t", quote = FALSE,
            row.names = FALSE)
