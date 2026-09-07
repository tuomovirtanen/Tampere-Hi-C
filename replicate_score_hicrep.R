###############################################################
## Pairwise HiCRep SCC calculation from .hic files
## GitHub version: https://github.com/TaoYang-dev/hicrep
###############################################################

## Install (run once)
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

if (!requireNamespace("HiCRep", quietly = TRUE))
    remotes::install_github("TaoYang-dev/hicrep")

library(hicrep)
library(tools)
#############################
## USER SETTINGS
#############################

hic_dir <- "../hic_files"

chromosome <- "1"      # chromosome to compare
resolution <- 100000   # 100 kb
normalization <- "KR"  # NONE, VC, VC_SQRT, KR

h <- 5                 # smoothing parameter
lbr <- 0               # lower bound (bp)
ubr <- 5000000         # upper bound (bp)

#############################
## FIND FILES
#############################

hic_files <- list.files(
    hic_dir,
    pattern = "\\.hic$",
    full.names = TRUE
)

sample_names <- file_path_sans_ext(basename(hic_files))

#############################
## RESULTS TABLE
#############################



#############################
## CACHE MATRICES
## (each .hic is read only once)
#############################

matrices <- vector("list", length(hic_files))

for(i in seq_along(hic_files)){

    message("Loading ", sample_names[i])

    matrices[[i]] <- hic2mat(
        file = hic_files[i],
        chromosome1 = chromosome,
        chromosome2 = chromosome,
        resol = resolution,
        method = normalization
    )
}

#############################
## PAIRWISE COMPARISONS
#############################
results <- data.frame(
    sample1 = character(),
    sample2 = character(),
    SCC = numeric(),
    stringsAsFactors = FALSE
)
for(i in 1:(length(hic_files)-1)){

    for(j in (i+1):length(hic_files)){

        message(sample_names[i], "  vs  ", sample_names[j])
        mat1 <- matrices[[i]]
        mat2 <- matrices[[j]]

        ## Crop both matrices to the smaller dimension
        common_size <- min(nrow(mat1), nrow(mat2))

        mat1 <- mat1[1:common_size, 1:common_size]
        mat2 <- mat2[1:common_size, 1:common_size]

        ## Remove bins that are entirely NA/NaN in either matrix
        bad1 <- apply(mat1, 1, function(x) sum(is.na(x) | is.nan(x)) > 50)
        bad2 <- apply(mat2, 1, function(x) sum(is.na(x) | is.nan(x)) > 50)

        remove_idx <- bad1 | bad2

        mat1 <- mat1[!remove_idx, !remove_idx, drop = FALSE]
        mat2 <- mat2[!remove_idx, !remove_idx, drop = FALSE]
        mat1[is.na(mat1)] <- 0
        mat2[is.na(mat2)] <- 0
        ## Skip if nothing remains
        if (nrow(mat1) < 2) {
            warning(sprintf("Skipping %s vs %s: too few bins after filtering.",
                            sample_names[i], sample_names[j]))
            next
        }
        scc <- get.scc(
            mat1,
            mat2,
            resol = resolution,
            h = h,
            lbr = lbr,
            ubr = ubr
        )

        results <- rbind(
            results,
            data.frame(
                sample1 = sample_names[i],
                sample2 = sample_names[j],
                SCC = scc$scc,
                stringsAsFactors = FALSE
            )
        )
    }
}

#############################
## SAVE RESULTS
#############################

write.table(
  results,
  file = "HiCRep_pairwise_SCC.tsv",
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

print(results)


library(dplyr)
library(tidyr)

# Read your data
df <- results

# Get all unique samples
samples <- sort(unique(c(df$sample1, df$sample2)))

# Create empty matrix
scc_matrix <- matrix(
  NA,
  nrow = length(samples),
  ncol = length(samples),
  dimnames = list(samples, samples)
)

# Diagonal = 1 (sample compared with itself)
diag(scc_matrix) <- 1

# Fill matrix
for (i in seq_len(nrow(df))) {
  s1 <- df$sample1[i]
  s2 <- df$sample2[i]
  scc <- df$SCC[i]

  scc_matrix[s1, s2] <- scc
  scc_matrix[s2, s1] <- scc
}

# Convert to data frame
scc_matrix <- as.data.frame(scc_matrix)

# Save matrix
write.table(
  scc_matrix,
  file = "SCC_matrix.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)