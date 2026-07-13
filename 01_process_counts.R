# 01_process_counts.R
# Purpose: Load raw counts and generate metadata matrix.

library(tidyverse)
library(here) # Best practice for reproducible file paths

# Define relative paths
raw_data_dir <- here("data", "raw_data")

# Reading the counts files
files <- list.files(path = raw_data_dir, pattern = "*.counts", full.names = TRUE)
sample_names <- tools::file_path_sans_ext(basename(files))

# Making the counts matrix
counts_list <- lapply(files, read.table, header=FALSE, row.names=1)
counts_matrix <- do.call(cbind, counts_list)
colnames(counts_matrix) <- sample_names

# Creation of metadata file
metadata_parts <- t(sapply(sample_names, function(s) {
  parts <- strsplit(s, "_")[[1]]
  if (length(parts) == 4) {
    return(parts)
  } else if (length(parts) == 3) {
    return(c(parts[1], parts[2], "NA", parts[3]))
  } else {
    stop(paste("Sample name error in ", s))
  }
}))

coldata <- data.frame(
  row.names = sample_names,
  genotype = metadata_parts[,1],
  timepoint = metadata_parts[,2],
  treatment = metadata_parts[,3],
  replicate = metadata_parts[,4],
  stringsAsFactors = FALSE
)

# Quality check
stopifnot(all(colnames(counts_matrix) == rownames(coldata)))

# Save processed inputs for the next script
dir.create(here("data", "processed"), recursive = TRUE, showWarnings = FALSE)
saveRDS(counts_matrix, here("data", "processed", "counts_matrix.rds"))
saveRDS(coldata, here("data", "processed", "coldata.rds"))