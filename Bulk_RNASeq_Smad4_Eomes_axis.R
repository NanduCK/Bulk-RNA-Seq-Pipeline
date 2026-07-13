#Bulk RNA-Seq script

#Loading libraries

library(DESeq2)
library(tidyverse)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)
library(pheatmap)
library(GEOquery)
library(oligo)
library(mogene10sttranscriptcluster.db)
library(limma)
library(fgsea)

setwd("D:/From_downloads/SMAD_KO_Exp/raw data/")

#reading the counts files
files <- list.files(pattern = "*.counts",full.names = TRUE)
sample_names <- tools::file_path_sans_ext(basename(files))


#Making the counts matrix
counts_list <- lapply(files, read.table, header=FALSE, row.names=1)

counts_matrix <- do.call(cbind,counts_list)
colnames(counts_matrix) <- sample_names


#CREATION of metadata file
metadata_parts <- t(sapply(sample_names, function(s) {
  parts <- strsplit(s,"_")[[1]]
  #we use strsplit to split the sample name at the _ , s is the argument of the function and [[1]] means that we are getting the vector alone instead of it as a list
  if(length(parts)==4) {
    #i.e genotype timepoint treatment replicate  
    return(parts)
  }
  else if(length(parts)==3){
    #i.e the treatment is N.A
    return(c(parts[1],parts[2],"NA", parts[3]))
  }
  else {
    stop(paste("Sample name error in ",s))
  }
  
}  ))

coldata <- data.frame(
  row.names = sample_names,
  genotype = metadata_parts[,1],
  timepoint = metadata_parts[,2],
  treatment = metadata_parts[,3],
  replicate = metadata_parts[,4],
  stringsAsFactors = FALSE
  
)

#checking to make sure row names match the count matrix, we should get TRUE if everything is correct
all(colnames(counts_matrix)==rownames(coldata))



####################################################################################################################################################

#DESeq


#Since genotype,timepoint and treatment are collinear(ie treatment is NA for 0d and 0hr,it creates a not fully ranked matrix).Therefore we collapse all the three variables into one
coldata$group <- factor(paste(coldata$genotype,coldata$timepoint,coldata$treatment, sep = "_"))

#Deseqdataset created using the new collapsed variable
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = coldata,
                              design = ~group)
#running DESeq
dds <- DESeq(dds)

#testing contrast
res1 <- results(dds, contrast = c("group", "Ctrl_0d_NA","Ctrl_0hr_NA" )) %>% as.data.frame() %>% na.omit()
res2 <- results(dds, contrast = c("group", "Ctrl_0d_NA","Ctrl_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res3 <- results(dds, contrast = c("group", "Ctrl_0d_NA","Ctrl_2hr_NoT" )) %>% as.data.frame() %>% na.omit()
res4 <- results(dds, contrast = c("group", "Ctrl_0d_NA","Ctrl_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res5 <- results(dds, contrast = c("group", "Ctrl_0d_NA","S34_0d_NA" )) %>% as.data.frame() %>% na.omit()
res6 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","Ctrl_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res7 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","Ctrl_2hr_NoT" )) %>% as.data.frame() %>% na.omit()
res8 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","Ctrl_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res9 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","EKO_0hr_NA" )) %>% as.data.frame() %>% na.omit()
res10 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","EKO_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res11 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","S34_0hr_NA" )) %>% as.data.frame() %>% na.omit()
res12 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","S34_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res13 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","S34_2hr_NoT" )) %>% as.data.frame() %>% na.omit()
res14 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","S34_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res15 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","S4_0hr_NA" )) %>% as.data.frame() %>% na.omit()
res16 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","S4_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res17 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","S4_2hr_NoT" )) %>% as.data.frame() %>% na.omit()
res18 <- results(dds, contrast = c("group", "Ctrl_0hr_NA","S4_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res19 <- results(dds, contrast = c("group", "Ctrl_12hr_NoT","Ctrl_12hr_T" )) %>% as.data.frame() %>% na.omit()
res20 <- results(dds, contrast = c("group", "Ctrl_12hr_NoT","Ctrl_2hr_NoT" )) %>% as.data.frame() %>% na.omit()
res21 <- results(dds, contrast = c("group", "Ctrl_12hr_NoT","Ctrl_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res22 <- results(dds, contrast = c("group", "Ctrl_12hr_NoT","EKO_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res23 <- results(dds, contrast = c("group", "Ctrl_12hr_NoT","S34_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res24 <- results(dds, contrast = c("group", "Ctrl_12hr_NoT","S4_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res25 <- results(dds, contrast = c("group", "Ctrl_12hr_NoT","S4_12hr_T" )) %>% as.data.frame() %>% na.omit()
res26 <- results(dds, contrast = c("group", "Ctrl_12hr_T","Ctrl_2hr_T" )) %>% as.data.frame() %>% na.omit()
res27 <- results(dds, contrast = c("group", "Ctrl_12hr_T","Ctrl_6hr_T" )) %>% as.data.frame() %>% na.omit()
res28 <- results(dds, contrast = c("group", "Ctrl_12hr_T","S34_12hr_T" )) %>% as.data.frame() %>% na.omit()
res29 <- results(dds, contrast = c("group", "Ctrl_12hr_T","S4_12hr_T" )) %>% as.data.frame() %>% na.omit()
res30 <- results(dds, contrast = c("group", "Ctrl_2hr_NoT","Ctrl_2hr_T" )) %>% as.data.frame() %>% na.omit()
res31 <- results(dds, contrast = c("group", "Ctrl_2hr_NoT","Ctrl_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res32 <- results(dds, contrast = c("group", "Ctrl_2hr_NoT","S34_2hr_NoT" )) %>% as.data.frame() %>% na.omit()
res33 <- results(dds, contrast = c("group", "Ctrl_2hr_T","Ctrl_6hr_T" )) %>% as.data.frame() %>% na.omit()
res34 <- results(dds, contrast = c("group", "Ctrl_2hr_T","S34_2hr_T" )) %>% as.data.frame() %>% na.omit()
res35 <- results(dds, contrast = c("group", "Ctrl_2hr_T","S4_2hr_T" )) %>% as.data.frame() %>% na.omit()
res36 <- results(dds, contrast = c("group", "Ctrl_6hr_NoT","Ctrl_6hr_T" )) %>% as.data.frame() %>% na.omit()
res37 <- results(dds, contrast = c("group", "Ctrl_6hr_NoT","EKO_6hr_T" )) %>% as.data.frame() %>% na.omit()
res38 <- results(dds, contrast = c("group", "Ctrl_6hr_NoT","S34_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res39 <- results(dds, contrast = c("group", "Ctrl_6hr_NoT","S4_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res40 <- results(dds, contrast = c("group", "Ctrl_6hr_T","EKO_6hr_T" )) %>% as.data.frame() %>% na.omit()
res41 <- results(dds, contrast = c("group", "Ctrl_6hr_T","S34_6hr_T" )) %>% as.data.frame() %>% na.omit()
res42 <- results(dds, contrast = c("group", "Ctrl_6hr_T","S4_6hr_T" )) %>% as.data.frame() %>% na.omit()
res43 <- results(dds, contrast = c("group", "EKO_0hr_NA","EKO_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res44 <- results(dds, contrast = c("group", "EKO_0hr_NA","EKO_6hr_T" )) %>% as.data.frame() %>% na.omit()
res45 <- results(dds, contrast = c("group", "EKO_0hr_NA","S34_0hr_NA" )) %>% as.data.frame() %>% na.omit()
res46 <- results(dds, contrast = c("group", "EKO_0hr_NA","S4_0hr_NA" )) %>% as.data.frame() %>% na.omit()
res47 <- results(dds, contrast = c("group", "EKO_12hr_NoT","S34_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res48 <- results(dds, contrast = c("group", "EKO_12hr_NoT","S4_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res49 <- results(dds, contrast = c("group", "EKO_6hr_T","S34_6hr_T" )) %>% as.data.frame() %>% na.omit()
res50 <- results(dds, contrast = c("group", "EKO_6hr_T","S4_6hr_T" )) %>% as.data.frame() %>% na.omit()
res51 <- results(dds, contrast = c("group", "S34_0d_NA","S34_0hr_NA" )) %>% as.data.frame() %>% na.omit()
res52 <- results(dds, contrast = c("group", "S34_0hr_NA","S34_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res53 <- results(dds, contrast = c("group", "S34_0hr_NA","S34_2hr_NoT" )) %>% as.data.frame() %>% na.omit()
res54 <- results(dds, contrast = c("group", "S34_0hr_NA","S4_0hr_NA" )) %>% as.data.frame() %>% na.omit()
res55 <- results(dds, contrast = c("group", "S34_0hr_NA","S34_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res56 <- results(dds, contrast = c("group", "S34_12hr_NoT","S34_12hr_T" )) %>% as.data.frame() %>% na.omit()
res57 <- results(dds, contrast = c("group", "S34_12hr_NoT","S34_2hr_NoT" )) %>% as.data.frame() %>% na.omit()
res58 <- results(dds, contrast = c("group", "S34_12hr_NoT","S34_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res59 <- results(dds, contrast = c("group", "S34_12hr_NoT","S4_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res60 <- results(dds, contrast = c("group", "S34_12hr_NoT","S4_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res61 <- results(dds, contrast = c("group", "S34_12hr_T","S34_2hr_T" )) %>% as.data.frame() %>% na.omit()
res62 <- results(dds, contrast = c("group", "S34_12hr_T","S34_6hr_T" )) %>% as.data.frame() %>% na.omit()
res63 <- results(dds, contrast = c("group", "S34_12hr_T","S4_12hr_T" )) %>% as.data.frame() %>% na.omit()
res64 <- results(dds, contrast = c("group", "S34_2hr_NoT","S34_2hr_T" )) %>% as.data.frame() %>% na.omit()
res65 <- results(dds, contrast = c("group", "S34_2hr_NoT","S34_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res66 <- results(dds, contrast = c("group", "S34_2hr_NoT","S4_2hr_NoT" )) %>% as.data.frame() %>% na.omit()
res67 <- results(dds, contrast = c("group", "S34_2hr_T","S34_6hr_T" )) %>% as.data.frame() %>% na.omit()
res68 <- results(dds, contrast = c("group", "S34_2hr_T","S4_2hr_T" )) %>% as.data.frame() %>% na.omit()
res69 <- results(dds, contrast = c("group", "S34_6hr_NoT","S34_6hr_T" )) %>% as.data.frame() %>% na.omit()
res70 <- results(dds, contrast = c("group", "S34_6hr_NoT","S4_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res71 <- results(dds, contrast = c("group", "S34_6hr_T","S4_6hr_T" )) %>% as.data.frame() %>% na.omit()
res72 <- results(dds, contrast = c("group", "S4_0hr_NA","S4_12hr_NoT" )) %>% as.data.frame() %>% na.omit()
res73 <- results(dds, contrast = c("group", "S4_0hr_NA","S4_2hr_NoT" )) %>% as.data.frame() %>% na.omit()
res74 <- results(dds, contrast = c("group", "S4_0hr_NA","S4_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res75 <- results(dds, contrast = c("group", "S4_12hr_NoT","S4_12hr_T" )) %>% as.data.frame() %>% na.omit()
res76 <- results(dds, contrast = c("group", "S4_12hr_NoT","S4_2hr_NoT" )) %>% as.data.frame() %>% na.omit()
res77 <- results(dds, contrast = c("group", "S4_12hr_NoT","S4_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res78 <- results(dds, contrast = c("group", "S4_12hr_T","S4_2hr_T" )) %>% as.data.frame() %>% na.omit()
res79 <- results(dds, contrast = c("group", "S4_12hr_T","S4_6hr_T" )) %>% as.data.frame() %>% na.omit()
res80 <- results(dds, contrast = c("group", "S4_2hr_NoT","S4_2hr_T" )) %>% as.data.frame() %>% na.omit()
res81 <- results(dds, contrast = c("group", "S4_2hr_NoT","S4_6hr_NoT" )) %>% as.data.frame() %>% na.omit()
res82 <- results(dds, contrast = c("group", "S4_2hr_T","S4_6hr_T" )) %>% as.data.frame() %>% na.omit()
res83 <- results(dds, contrast = c("group", "Ctrl_2hr_NoT","S4_2hr_NoT" )) %>% as.data.frame() %>% na.omit()



#Conversion of MGI Ids to Entrez ids to create ranked lists

#Create an empty list to hold the processed data frames
processed_res_list <- list()

#Loop from 1 to 83
for (i in 1:83) {
  #getting the name of the dataframe dynamically
  df_name <- paste0("res",i)
  #checking if such a dataframe exists
  if(exists(df_name) && is.object(get(df_name))){
    message(paste("Processing:",df_name))
    #getting the dataframes and then getting the entrezids
    current_df <-as.data.frame(get(df_name)) %>%
      tibble::rownames_to_column(var = "original_rownames")
    current_df$geneid <- sub(".*GeneID:(\\d+).*", "\\1", current_df$original_rownames)
    #storing it all into the list
    #[[...]] this is called list accessor that lets you use a variable to specify the name of the element you want to create or access.
    processed_res_list[[df_name]] <- current_df
  }else{
    message(paste("Skipping",df_name," (does not exist)"))
  }
  
}
#doing the actual mapping
mapped_res_list<- list()

for (contrast_name in names(processed_res_list)) {
  temp_df <- processed_res_list[[contrast_name]]
  temp_df$SYMBOL <- mapIds(org.Mm.eg.db,
                           keys = temp_df$geneid,
                           column = "SYMBOL",
                           keytype = "ENTREZID",
                           multiVals = "list")
  mapped_res_list[[contrast_name]] <- temp_df
  
  message(paste("Mapped the entrez ids to symbols successfully for", contrast_name))
}

#creating ranked list

output_dir <- "D:/From_downloads/SMAD_KO_Exp/results/GSEA_input"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  message(paste("Created output directory:", output_dir))
}

stat_col_name<- "stat"

for (contrast_name in names(mapped_res_list)) {
  temp_df <- mapped_res_list[[contrast_name]]
  ranked_df<-temp_df %>% 
    #removing genes not mapped to a symbol
    filter(!is.na(SYMBOL)) %>%
    group_by(SYMBOL) %>%
    #For each symbol, keep only the row with the largest *absolute* stat value.
    slice_max(order_by = abs(.data[[stat_col_name]]), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(SYMBOL, all_of(stat_col_name)) %>%
    #SYMBOL should be character and stat should be numeric so flattening them using mutate
    mutate(
      SYMBOL = as.character(SYMBOL),
      !!sym(stat_col_name) := as.numeric(.data[[stat_col_name]])
    ) %>%
    arrange(desc(.data[[stat_col_name]]))
  
  output_filename <- file.path(output_dir, paste0(contrast_name, ".rnk"))
  
  write.table(ranked_df, 
              file = output_filename,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  message(paste("Saved:", output_filename))
  
}



#getting normalised counts matrix
size_factors <- estimateSizeFactors(dds) %>% as.data.frame()

norm_counts <- counts(dds, normalized= TRUE) %>%as.data.frame()

#COnverting MGI ids to entrez ids to map 
norm_counts$geneid <- row.names(norm_counts)
norm_counts <- norm_counts[,c(90,1:89)]
norm_counts$geneid <- sub(".*GeneID:(\\d+).*", "\\1", norm_counts$geneid)


###################################################################################################################################################

#for gse39152 - Wakim et al(2012) - Brain TRM GSEA

# 1. Define the gene lists- Taken directly from the paper
brain_up <- c("Icos", "Skil", "Hsph1", "Ccl3", "Pdcd1", "Gadd45b", "Il12rb2", "Dtx4", "Rgs1", 
              "Dnajb1", "Nfil3", "Dusp4", "Fosl2", "Pmepa1", "Ppp1r15a", "Coro2a", "Tigit", 
              "Egr1", "Bag3", "Rgs2", "Phlda1", "Sik1", "Vps37b", "Il21r", "Il4ra", "Tnfrsf1b", 
              "Dusp2", "Slc3a2", "Neurl3", "Mxd1", "Adam19", "Pik3ap1", "Nr4a3", "Cxcl10", 
              "Gm10008", "Ifi44", "Hspa1a", "Rhob")

brain_down <- c("A430078G23Rik", "Dnahc8", "Arl5c", "Ces2", "Cdc14b", "Wfkkmn2", "Bin2", 
                "Rasa3", "Klhl6", "Arl4c", "Elovl7", "Bclp2")

# 2. Format as a GMT line: Name \t Description \t Genes...
line1 <- paste(c("BRAIN_TRM_UP", "Genes_upregulated_in_Brain_CD103plus", brain_up), collapse = "\t")
line2 <- paste(c("BRAIN_TRM_DN", "Genes_downregulated_in_Brain_CD103plus", brain_down), collapse = "\t")

# 3. Write to file
writeLines(c(line1, line2), "Wakim_Brain_TRM.gmt")






#for gse70183 - Mackay et al(2016)

Mackaycounts <- read.table("GSE70813/GSE70813_Supp_Raw_Counts.txt.gz", header=TRUE, row.names=1)


Mackaycounts_clean <- Mackaycounts %>%
  # 1. Remove rows where Symbol is NA or just an empty string
  filter(!is.na(Symbol), Symbol != "", Symbol != " ") %>% 
  
  # 2. Proceed with your deduplication logic
  mutate(TotalExp = rowSums(dplyr::select(., where(is.numeric)))) %>% 
  group_by(Symbol) %>%
  slice_max(order_by = TotalExp, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  
  # 3. Clean up and finalize
  dplyr::select(-TotalExp) %>%
  column_to_rownames("Symbol")



Mtarget_columns <- c(
  "C2_C5N42ANXX_CGATGTAT.bam", "D2_C5N42ANXX_GGCTACAT.bam", "Spleen_TCM_Sort1_C6GUTANXX_AGTCAACA_L003.bam" , "Spleen_TCM_Sort2_C6GUTANXX_CCGTCCCG_L003.bam", # Spleen TCM
  "C3_C5N42ANXX_TTAGGCAT.bam", "D3_C5N42ANXX_CTTGTAAT.bam", "Spleen_TEM_Sort1_C6GUTANXX_AGTTCCGT_L003.bam", "Spleen_TEM_Sort2_C6GUTANXX_GTCCGCAC_L003.bam", # Spleen TEM
  "C4_C5N42ANXX_TGACCAAT.bam", "D4_C5N42ANXX_AGTCAACA.bam", # Liver TEM
  "C5_C5N42ANXX_ACAGTGAT.bam", "D5_C5N42ANXX_AGTTCCGT.bam" #Liver TRM
)

Mcol_data <- data.frame(
  row.names = Mtarget_columns,
  condition = factor(c(
    rep("Spleen_TCM", 4), 
    rep("Spleen_TEM", 4), 
    rep("Liver_TEM", 2), 
    rep("Liver_TRM", 2)
  ))
)

Mackaydds <- DESeqDataSetFromMatrix(countData = Mackaycounts_clean[,Mtarget_columns],
                                    colData = Mcol_data,
                                    design = ~ condition)
Mackaydds <- DESeq(Mackaydds)

Mres1 <- results(Mackaydds, contrast=c("condition", "Liver_TRM", "Liver_TEM"))  %>% as.data.frame() %>% na.omit()
Mres2 <- results(Mackaydds, contrast=c("condition", "Liver_TEM", "Spleen_TEM"))  %>% as.data.frame() %>% na.omit()
Mres3 <- results(Mackaydds, contrast=c("condition", "Spleen_TEM", "Spleen_TCM"))  %>% as.data.frame() %>% na.omit()
Mres4 <- results(Mackaydds, contrast = c("condition", "Spleen_TCM", "Spleen_TEM")) %>% as.data.frame() %>% na.omit()


# Extract only the UPREGULATED genes (LFC > 1 and padj < 0.05)
get_up_genes <- function(res) {
  res_df <- as.data.frame(res)
  rownames(res_df[which(res_df$log2FoldChange > 1 & res_df$padj < 0.05), ])
}

LiverTRM_genes <- get_up_genes(Mres1)
LiverTEM_genes <- get_up_genes(Mres2)
SpleenTEM_genes <- get_up_genes(Mres3)
SpleenTCM_genes <- get_up_genes(Mres4)

gmt_lines <- c(
  paste(c("MACKAY_LIVER_TRM_UP", "Upregulated_in_Liver_TRM_vs_TEM", LiverTRM_genes), collapse = "\t"),
  paste(c("MACKAY_LIVER_TEM_UP", "Upregulated_in_Liver_TEM_vs_Spleen_TEM", LiverTEM_genes), collapse = "\t"),
  paste(c("MACKAY_SPLEEN_TEM_UP", "Upregulated_in_Spleen_TEM_vs_Spleen_TCM", SpleenTEM_genes), collapse = "\t"),
  paste(c("MACKAY_SPLEEN_TCM_UP", "Upregulated_in_Spleen_TCM_vs_Spleen_TEM", SpleenTCM_genes), collapse = "\t")
)

writeLines(gmt_lines, "Mackay_Four_Populations.gmt")


########################################################################################################################################


#Plotting

res19_vector <- setNames(res19$stat, rownames(res19))
res22_vector <- setNames(res22$stat, rownames(res22))
res23_vector <- setNames(res23$stat, rownames(res23))
res24_vector <- setNames(res24$stat, rownames(res24))


multi_sample_list <- list(
  "Ctrl_12hr_NoT vs Ctrl_12hr_T" = res19_vector,
  "Ctrl_12hr_NoT vs EKO_12hr_NoT" = res22_vector,
  "Ctrl_12hr_NoT vs S34_12hr_NoT" = res23_vector,
  "Ctrl_12hr_NoT vs S4_12hr_NoT" = res24_vector
)

#Creating a combined graph of mackay and wakim

combined_pathways <- list(
  "Wakim_Brain_TRM" = brain_up, 
  "Mackay_Liver_TRM" = LiverTRM_genes,
  "Mackay_Spleen_TCM" = SpleenTCM_genes,
  "Mackay_Spleen_TEM" = SpleenTEM_genes
)

  
  
fgsea_res19 <- fgsea(pathways = combined_pathways, 
                       stats = res19_vector, 
                       minSize = 5, 
                       maxSize = 500)


#optional(
# Check overlap for each pathway
overlap_summary <- lapply(names(combined_pathways), function(name) {
  genes_in_set <- combined_pathways[[name]]
  genes_in_ranks <- names(res19_vector)
  matches <- intersect(genes_in_set, genes_in_ranks)
  
  data.frame(Pathway = name, 
             Set_Size = length(genes_in_set), 
             Matches_Found = length(matches))
})

do.call(rbind, overlap_summary)

#)

# We use 'plotEnrichmentData' to get the raw points for the curve

# This creates a data frame of the running ES for each rank
plot_list_res19 <- lapply(names(combined_pathways), function(pw_name) {
  pathway_genes <- combined_pathways[[pw_name]]
  common_genes <- intersect(pathway_genes, names(res19_vector))
  
  if (length(common_genes) < 2) {
    message(paste("Skipping", pw_name, "- only", length(common_genes), "matches found."))
    return(NULL)
  }
  gsea_res <- tryCatch({
    plotEnrichmentData(pathway_genes, res19_vector)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(gsea_res) || !is.list(gsea_res)) {
    return(NULL)
  }
  curve_df <- data.frame(
    rank = gsea_res$curve$rank,
    ES   = gsea_res$curve$ES,
    Pathway = pw_name,
    stringsAsFactors = FALSE
  )
  ticks_df <- data.frame(
    rank = gsea_res$ticks,
    Pathway = pw_name,
    stringsAsFactors = FALSE
  )
  return(list(curve = curve_df, ticks = ticks_df))
})
  
# Combine all curves into one data frame for ggplot
plot_list_res19 <- compact(plot_list_res19)
all_curves_res19 <- bind_rows(map(plot_list_res19, "curve"))
all_ticks_res19  <- bind_rows(map(plot_list_res19, "ticks"))
all_ticks_res19 <- all_ticks_res19 %>%
  rename(rank = rank.rank) # Change 'rank.rank' back to 'rank'



all_ticks_res19 <- all_ticks_res19 %>%
  mutate(y_start = -0.1 - (as.numeric(as.factor(Pathway)) * 0.05),
         y_end   = y_start - 0.04)


# 1. Check if we actually have data
print(nrow(all_ticks_res19)) 

# 2. Check the column names
print(colnames(all_ticks_res19))



ggplot() +
  # Layer 1: The Enrichment Curves
  geom_line(data = all_curves_res19, aes(x = rank, y = ES, color = Pathway), size = 1.2) +
  
  # Layer 2: The Barcode Ticks (now using the columns we explicitly built)
  geom_segment(data = all_ticks_res19, 
               aes(x = rank, xend = rank, y = y_start, yend = y_end, color = Pathway)) +
  
  # Layer 3: The Zero Line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Formatting
  scale_color_manual(values = c(
    "Wakim_Brain_TRM"   = "#1f78b4", # Blue
    "Mackay_Liver_TRM"  = "#e31a1c", # Red
    "Mackay_Spleen_TEM" = "#33a02c", # Green
    "Mackay_Spleen_TCM" = "#ff7f00"  # Orange
  )) +
  theme_minimal() +
  labs(title = "Ctrl_12hr_NoT vsCtrl_12hr_T",
       subtitle = "Comparing Brain (Wakim) vs. Liver/Spleen (Mackay)",
       x = "Rank in Gene List",
       y = "Enrichment Score (ES)") +
  theme(legend.position = "bottom")










#Multi Sample Plot

res24$geneid <- rownames(res24) 
res24$geneid <- sub(".*GeneID:(\\d+).*", "\\1", res24$geneid)
res24$SYMBOL <- mapIds(org.Mm.eg.db,
                       keys = res24$geneid,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "list")

res24 <- res24[order(-res24$stat),] %>% filter(!is.na(SYMBOL)) %>%group_by(SYMBOL) %>%
  slice_max(order_by = abs(.data[[stat_col_name]]), n = 1, with_ties = FALSE) %>%
  ungroup() 
rownames(res24) <- res24$SYMBOL

shared_data <- read.csv("D:/From_downloads/shared_genes.csv")

# Create Smad4 Ranking (using LFC and p-value)
# We handle p-values of 0 by adding a tiny offset to avoid infinity
res_smad4_vector <- shared_data %>%
  mutate(stat = sign(lfc_S4) * -log10(pval_S4 + 1e-300)) %>%
  filter(!is.na(stat), !is.na(SYMBOL)) %>%
  { setNames(.$stat, .$SYMBOL) }

# Create Eomes Ranking
res_eomes_vector <- shared_data %>%
  mutate(stat = sign(lfc_Eomes) * -log10(pval_Eomes + 1e-300)) %>%
  filter(!is.na(stat), !is.na(SYMBOL)) %>%
  { setNames(.$stat, .$SYMBOL) }



#faceted plot
multi_sample_curves <- map_df(names(combined_pathways), function(pw_name) {
  
  map_df(names(multi_sample_list), function(samp_name) {
    
    # Get the specific pathway and sample vector
    pathway_genes <- combined_pathways[[pw_name]]
    current_ranks <- multi_sample_list[[samp_name]]
    
    # Run the curve calculation
    gsea_res <- tryCatch({
      plotEnrichmentData(pathway_genes, current_ranks)
    }, error = function(e) return(NULL))
    
    if (is.null(gsea_res)) return(NULL)
    
    # Return cleaned curve data
    tibble(
      rank = as.numeric(gsea_res$curve$rank),
      ES   = as.numeric(gsea_res$curve$ES),
      Pathway = pw_name,
      Sample  = samp_name
    )
  })
})

ggplot(multi_sample_curves, aes(x = rank, y = ES, color = Sample)) +
  # Layer 1: The Curves (No Ticks as requested)
  geom_line(size = 1) +
  
  # Layer 2: The Zero Line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Layer 3: Faceting by Gene Set
  facet_wrap(~Pathway, scales = "free_y") +
  
  # Formatting
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "grey80", fill = NA),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Multi-Sample GSEA Comparison",
    x = "Rank in Gene List",
    y = "Enrichment Score (ES)",
    color = "Sample ID"
  )





#Seperate plots

pathway_names <- unique(multi_sample_curves$Pathway)
gsea_plot_list <- list()

for (pw in pathway_names) {
  
  # Filter data for just this one pathway
  plot_data_sub <- multi_sample_curves %>% filter(Pathway == pw)
  
  # Create the plot for this specific gene set
  p <- ggplot(plot_data_sub, aes(x = rank, y = ES, color = Sample)) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    theme_minimal() +
    labs(
      title = paste("GSEA Profile:", pw),
      subtitle = "Overlaid Samples Comparison",
      x = "Rank in Gene List",
      y = "Enrichment Score (ES)"
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14)
    )
  
  # Store the plot in the list
  gsea_plot_list[[pw]] <- p
}

gsea_plot_list[["Mackay_Liver_TRM"]]
gsea_plot_list[["Wakim_Brain_TRM"]]
gsea_plot_list[["Mackay_Spleen_TCM"]]
gsea_plot_list[["Mackay_Spleen_TEM"]]

#smad and eomes genes

smad4_eomes_list <- list(
  "Smad4_Controlled" = res_smad4_vector,
  "Eomes_Controlled" = res_eomes_vector
)

smad_eomes_curves <- map_df(names(combined_pathways), function(pw_name) {
  map_df(names(smad4_eomes_list), function(samp_name) {
    
    pathway_genes <- combined_pathways[[pw_name]]
    current_ranks <- smad4_eomes_list[[samp_name]]
    
    # Run the curve calculation with error handling
    gsea_res <- tryCatch({
      plotEnrichmentData(pathway_genes, current_ranks)
    }, error = function(e) return(NULL))
    
    if (is.null(gsea_res) || is.null(gsea_res$curve)) return(NULL)
    
    # Return cleaned curve data with explicit names to avoid 'rank.rank'
    tibble(
      rank = as.numeric(gsea_res$curve$rank),
      ES   = as.numeric(gsea_res$curve$ES),
      Pathway = pw_name,
      Factor  = samp_name
    )
  })
})

shared_plots <- list()

for (pw in names(combined_pathways)) {
  
  plot_data_sub <- smad_eomes_curves %>% filter(Pathway == pw)
  
  if (nrow(plot_data_sub) == 0) next
  
  q <- ggplot(plot_data_sub, aes(x = rank, y = ES, color = Factor)) +
    geom_line(size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Smad4_Controlled" = "#D55E00", "Eomes_Controlled" = "#0072B2")) +
    theme_minimal() +
    labs(
      title = paste("Enrichment Profile:", pw),
      subtitle = "Overlay: Smad4 vs. Eomes Ranking",
      x = "Rank (Shared Genes)",
      y = "Enrichment Score (ES)"
    ) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  
  shared_plots[[pw]] <- q
}

shared_plots[["Mackay_Liver_TRM"]]
shared_plots[["Wakim_Brain_TRM"]]
shared_plots[["Mackay_Spleen_TCM"]]
shared_plots[["Mackay_Spleen_TEM"]]


for (name in names(gsea_plot_list)) {
  file_name <- paste0("GSEA_Plot_", name, ".png")
  ggsave(file_name, plot = gsea_plot_list[[name]], width = 10, height = 6, dpi = 300)
}

for (name in names(shared_plots)) {
  file_name <- paste0("GSEA_Plot_shared_", name, ".png")
  ggsave(file_name, plot = shared_plots[[name]], width = 10, height = 6, dpi = 300)
}


gsea_stats_report <- map_df(names(multi_sample_list), function(samp_name) {
  
  res <- fgsea(pathways = combined_pathways, 
               stats = multi_sample_list[[samp_name]],
               minSize = 10,
               maxSize = 500)
  
  # Add the sample name so we know which row is which
  res %>% 
    mutate(Sample = samp_name) %>%
    as_tibble()
})

final_report <- gsea_stats_report %>%
  dplyr::select(
    Sample,
    Pathway = pathway,
    Size = size,
    ES = ES,
    NES = NES,
    P_Value = pval,
    FDR_q_val = padj,  # 'padj' in fgsea is the FDR (BH method)
    leadingEdge
  ) %>%
  arrange(Sample, desc(NES))

write.csv(final_report %>% dplyr::select(-leadingEdge), "GSEA_MackayandWakim_Summary.csv", row.names = FALSE)




















