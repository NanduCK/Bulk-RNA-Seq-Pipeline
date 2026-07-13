# 02_run_deseq2.R
# Purpose: Execute differential expression analysis and generate ranked lists for GSEA.

library(DESeq2)
library(tidyverse)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(here)

counts_matrix <- readRDS(here("data", "processed", "counts_matrix.rds"))
coldata <- readRDS(here("data", "processed", "coldata.rds"))

# Collapse variables for design formula
coldata$group <- factor(paste(coldata$genotype, coldata$timepoint, coldata$treatment, sep = "_"))

dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = coldata, design = ~group)
dds <- DESeq(dds)

# Define specific contrast pairs 
contrasts_to_test <- list(
  c("Ctrl_0d_NA", "Ctrl_0hr_NA"),
  c("Ctrl_0d_NA", "Ctrl_12hr_NoT"),
  c("Ctrl_0d_NA", "Ctrl_2hr_NoT"),
  c("Ctrl_0d_NA", "Ctrl_6hr_NoT"),
  c("Ctrl_0d_NA", "S34_0d_NA"),
  c("Ctrl_0hr_NA", "Ctrl_12hr_NoT"),
  c("Ctrl_0hr_NA", "Ctrl_2hr_NoT"),
  c("Ctrl_0hr_NA", "Ctrl_6hr_NoT"),
  c("Ctrl_0hr_NA", "EKO_0hr_NA"),
  c("Ctrl_0hr_NA", "EKO_12hr_NoT"),
  c("Ctrl_0hr_NA", "S34_0hr_NA"),
  c("Ctrl_0hr_NA", "S34_12hr_NoT"),
  c("Ctrl_0hr_NA", "S34_2hr_NoT"),
  c("Ctrl_0hr_NA", "S34_6hr_NoT"),
  c("Ctrl_0hr_NA", "S4_0hr_NA"),
  c("Ctrl_0hr_NA", "S4_12hr_NoT"),
  c("Ctrl_0hr_NA", "S4_2hr_NoT"),
  c("Ctrl_0hr_NA", "S4_6hr_NoT"),
  c("Ctrl_12hr_NoT", "Ctrl_12hr_T"),
  c("Ctrl_12hr_NoT", "Ctrl_2hr_NoT"),
  c("Ctrl_12hr_NoT", "Ctrl_6hr_NoT"),
  c("Ctrl_12hr_NoT", "EKO_12hr_NoT"),
  c("Ctrl_12hr_NoT", "S34_12hr_NoT"),
  c("Ctrl_12hr_NoT", "S4_12hr_NoT"),
  c("Ctrl_12hr_NoT", "S4_12hr_T"),
  c("Ctrl_12hr_T", "Ctrl_2hr_T"),
  c("Ctrl_12hr_T", "Ctrl_6hr_T"),
  c("Ctrl_12hr_T", "S34_12hr_T"),
  c("Ctrl_12hr_T", "S4_12hr_T"),
  c("Ctrl_2hr_NoT", "Ctrl_2hr_T"),
  c("Ctrl_2hr_NoT", "Ctrl_6hr_NoT"),
  c("Ctrl_2hr_NoT", "S34_2hr_NoT"),
  c("Ctrl_2hr_T", "Ctrl_6hr_T"),
  c("Ctrl_2hr_T", "S34_2hr_T"),
  c("Ctrl_2hr_T", "S4_2hr_T"),
  c("Ctrl_6hr_NoT", "Ctrl_6hr_T"),
  c("Ctrl_6hr_NoT", "EKO_6hr_T"),
  c("Ctrl_6hr_NoT", "S34_6hr_NoT"),
  c("Ctrl_6hr_NoT", "S4_6hr_NoT"),
  c("Ctrl_6hr_T", "EKO_6hr_T"),
  c("Ctrl_6hr_T", "S34_6hr_T"),
  c("Ctrl_6hr_T", "S4_6hr_T"),
  c("EKO_0hr_NA", "EKO_12hr_NoT"),
  c("EKO_0hr_NA", "EKO_6hr_T"),
  c("EKO_0hr_NA", "S34_0hr_NA"),
  c("EKO_0hr_NA", "S4_0hr_NA"),
  c("EKO_12hr_NoT", "S34_12hr_NoT"),
  c("EKO_12hr_NoT", "S4_12hr_NoT"),
  c("EKO_6hr_T", "S34_6hr_T"),
  c("EKO_6hr_T", "S4_6hr_T"),
  c("S34_0d_NA", "S34_0hr_NA"),
  c("S34_0hr_NA", "S34_12hr_NoT"),
  c("S34_0hr_NA", "S34_2hr_NoT"),
  c("S34_0hr_NA", "S4_0hr_NA"),
  c("S34_0hr_NA", "S34_6hr_NoT"),
  c("S34_12hr_NoT", "S34_12hr_T"),
  c("S34_12hr_NoT", "S34_2hr_NoT"),
  c("S34_12hr_NoT", "S34_6hr_NoT"),
  c("S34_12hr_NoT", "S4_12hr_NoT"),
  c("S34_12hr_NoT", "S4_6hr_NoT"),
  c("S34_12hr_T", "S34_2hr_T"),
  c("S34_12hr_T", "S34_6hr_T"),
  c("S34_12hr_T", "S4_12hr_T"),
  c("S34_2hr_NoT", "S34_2hr_T"),
  c("S34_2hr_NoT", "S34_6hr_NoT"),
  c("S34_2hr_NoT", "S4_2hr_NoT"),
  c("S34_2hr_T", "S34_6hr_T"),
  c("S34_2hr_T", "S4_2hr_T"),
  c("S34_6hr_NoT", "S34_6hr_T"),
  c("S34_6hr_NoT", "S4_6hr_NoT"),
  c("S34_6hr_T", "S4_6hr_T"),
  c("S4_0hr_NA", "S4_12hr_NoT"),
  c("S4_0hr_NA", "S4_2hr_NoT"),
  c("S4_0hr_NA", "S4_6hr_NoT"),
  c("S4_12hr_NoT", "S4_12hr_T"),
  c("S4_12hr_NoT", "S4_2hr_NoT"),
  c("S4_12hr_NoT", "S4_6hr_NoT"),
  c("S4_12hr_T", "S4_2hr_T"),
  c("S4_12hr_T", "S4_6hr_T"),
  c("S4_2hr_NoT", "S4_2hr_T"),
  c("S4_2hr_NoT", "S4_6hr_NoT"),
  c("S4_2hr_T", "S4_6hr_T"),
  c("Ctrl_2hr_NoT", "S4_2hr_NoT")
)

# Run DESeq2 results over the list
deseq_results_list <- lapply(contrasts_to_test, function(pair) {
  results(dds, contrast = c("group", pair[1], pair[2])) %>%
    as.data.frame() %>%
    na.omit() %>%
    rownames_to_column(var = "original_rownames") %>%
    mutate(geneid = sub(".*GeneID:(\\d+).*", "\\1", original_rownames))
})

# Name the list elements logically
names(deseq_results_list) <- sapply(contrasts_to_test, function(x) paste(x[1], "vs", x[2], sep="_"))

# Mapping and Ranking
output_dir <- here("results", "GSEA_input")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stat_col_name <- "stat"
mapped_res_list <- list()

for (contrast_name in names(deseq_results_list)) {
  temp_df <- deseq_results_list[[contrast_name]]
  
  # Note: 'first' is preferred over 'list' here to keep the dataframe flat
  temp_df$SYMBOL <- mapIds(org.Mm.eg.db,
                           keys = temp_df$geneid,
                           column = "SYMBOL",
                           keytype = "ENTREZID",
                           multiVals = "first") 
  
  mapped_res_list[[contrast_name]] <- temp_df
  
  # Create ranked list
  ranked_df <- temp_df %>% 
    filter(!is.na(SYMBOL)) %>%
    group_by(SYMBOL) %>%
    slice_max(order_by = abs(.data[[stat_col_name]]), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(SYMBOL, all_of(stat_col_name)) %>%
    mutate(SYMBOL = as.character(SYMBOL),
           !!sym(stat_col_name) := as.numeric(.data[[stat_col_name]])) %>%
    arrange(desc(.data[[stat_col_name]]))
  
  output_filename <- file.path(output_dir, paste0(contrast_name, ".rnk"))
  write.table(ranked_df, file = output_filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# Getting normalized counts matrix safely without magic numbers
size_factors <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "geneid") %>%
  dplyr::select(geneid, everything()) %>% # Safely pulls geneid to the front
  mutate(geneid = sub(".*GeneID:(\\d+).*", "\\1", geneid))

saveRDS(mapped_res_list, here("data", "processed", "mapped_res_list.rds"))