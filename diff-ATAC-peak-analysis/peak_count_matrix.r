
# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript myscript.R <input_path>", call. = FALSE)
}

# Assign the first argument to 'input'
input <- args[1]

cat("Input path:", input, "\n")


library(data.table)
# 1. Get all peak count files
files <- list.files(input,pattern = "*peak_counts.bed.gz$")

# 2. Extract sample names (remove pattern suffix)
sample_names <- sub("_peak_counts\\.bed\\.gz$", "", files)

count_matrix <- matrix(NA, nrow = 0, ncol = length(sample_names) + 1)
# 3. Read all files into a named list
count_list <- lapply(paste0(input,"/",files), function(f) {
  dt <- fread(cmd = paste("zcat", f), select = c(4, 5), col.names = c("index", "count"))
  return(dt)
})
names(count_list) <- sample_names


count_matrix <- matrix(NA, nrow = nrow(count_list[[1]]), ncol = length(sample_names))
colnames(count_matrix) <- sample_names
for(s in sample_names){
    count_matrix[count_list[[s]]$index,s] <- count_list[[s]]$count
}

bed <- data.frame(fread(cmd = paste("zcat", paste0(input,"/",files[1])), select = c(1,2,3), col.names = c("chrom", "start","stop")))
# 8. Save to file if desired
write.table(cbind(bed,count_matrix), paste0(input,"/peak_count_matrix.tsv"), quote = FALSE, sep = "\t",row.names = FALSE)
