# 03_gsea_analysis.R
# Purpose: Integrate external datasets (GSE39152, GSE70183), perform GSEA, and generate plots.

library(tidyverse)
library(limma)
library(DESeq2)
library(fgsea)
library(oligo)
library(here)

# Create output directories for our results
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)

# Load the mapped DESeq2 results from Module 2
mapped_res_list <- readRDS(here("data", "processed", "mapped_res_list.rds"))

################################################################################
# 1. WAKIM DATASET (GSE39152) PROCESSING
################################################################################

wakim_exdir <- here("data", "external", "GSE39152_extracted_cel")
dir.create(wakim_exdir, recursive = TRUE, showWarnings = FALSE)

# Assuming tar is already extracted or handled. If not:
# untar(here("data", "external", "GSE39152_RAW.tar"), exdir = wakim_exdir)

cel_files_wakim <- list.files(wakim_exdir, pattern = ".CEL.gz", full.names = TRUE)
raw_data_wakim <- read.celfiles(cel_files_wakim)
norm_data_wakim <- rma(raw_data_wakim)
exp_matrix_wakim <- exprs(norm_data_wakim)

# Avoid magic numbers by extracting exactly 10 known columns (Brain TRM and Spleen)
target_wakim_samples <- colnames(exp_matrix_wakim)[4:13] 
exp_matrix_wakim_subset <- exp_matrix_wakim[, target_wakim_samples]

groups_wakim <- factor(c(rep("Spleen_Mem", 5), rep("Brain_TRM", 5)))
design_wakim <- model.matrix(~0 + groups_wakim)
colnames(design_wakim) <- c("Brain_TRM", "Spleen_Mem")

fit_wakim <- lmFit(exp_matrix_wakim_subset, design_wakim)
contrast.matrix_wakim <- makeContrasts(Brain_TRM - Spleen_Mem, levels=design_wakim)
fit_wakim2 <- eBayes(contrasts.fit(fit_wakim, contrast.matrix_wakim))

# Filter for UP-REGULATED genes in Brain TRM
brain_trm_results_wakim <- topTable(fit_wakim2, coef=1, number=Inf, adjust="fdr")
brain_trm_up <- brain_trm_results_wakim %>%
  filter(logFC > 1 & adj.P.Val < 0.05) %>%
  pull(SYMBOL)

# Hardcoded subsets from original script for GMTs
brain_up <- c("Icos", "Skil", "Hsph1", "Ccl3", "Pdcd1", "Gadd45b", "Il12rb2", "Dtx4", "Rgs1", 
              "Dnajb1", "Nfil3", "Dusp4", "Fosl2", "Pmepa1", "Ppp1r15a", "Coro2a", "Tigit", 
              "Egr1", "Bag3", "Rgs2", "Phlda1", "Sik1", "Vps37b", "Il21r", "Il4ra", "Tnfrsf1b", 
              "Dusp2", "Slc3a2", "Neurl3", "Mxd1", "Adam19", "Pik3ap1", "Nr4a3", "Cxcl10", 
              "Gm10008", "Ifi44", "Hspa1a", "Rhob")
brain_down <- c("A430078G23Rik", "Dnahc8", "Arl5c", "Ces2", "Cdc14b", "Wfkkmn2", "Bin2", 
                "Rasa3", "Klhl6", "Arl4c", "Elovl7", "Bclp2")

writeLines(c(
  paste(c("BRAIN_TRM_UP", "Genes_upregulated_in_Brain_CD103plus", brain_up), collapse = "\t"),
  paste(c("BRAIN_TRM_DN", "Genes_downregulated_in_Brain_CD103plus", brain_down), collapse = "\t")
), here("results", "tables", "Wakim_Brain_TRM.gmt"))

################################################################################
# 2. MACKAY DATASET (GSE70813) PROCESSING
################################################################################

Mackaycounts <- read.table(here("data", "external", "GSE70813_Supp_Raw_Counts.txt.gz"), header=TRUE, row.names=1)

Mackaycounts_clean <- Mackaycounts %>%
  filter(!is.na(Symbol), Symbol != "", Symbol != " ") %>% 
  mutate(TotalExp = rowSums(dplyr::select(., where(is.numeric)))) %>% 
  group_by(Symbol) %>%
  slice_max(order_by = TotalExp, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-TotalExp) %>%
  column_to_rownames("Symbol")

Mtarget_columns <- c(
  "C2_C5N42ANXX_CGATGTAT.bam", "D2_C5N42ANXX_GGCTACAT.bam", "Spleen_TCM_Sort1_C6GUTANXX_AGTCAACA_L003.bam", "Spleen_TCM_Sort2_C6GUTANXX_CCGTCCCG_L003.bam", 
  "C3_C5N42ANXX_TTAGGCAT.bam", "D3_C5N42ANXX_CTTGTAAT.bam", "Spleen_TEM_Sort1_C6GUTANXX_AGTTCCGT_L003.bam", "Spleen_TEM_Sort2_C6GUTANXX_GTCCGCAC_L003.bam", 
  "C4_C5N42ANXX_TGACCAAT.bam", "D4_C5N42ANXX_AGTCAACA.bam", 
  "C5_C5N42ANXX_ACAGTGAT.bam", "D5_C5N42ANXX_AGTTCCGT.bam" 
)

Mcol_data <- data.frame(
  row.names = Mtarget_columns,
  condition = factor(c(rep("Spleen_TCM", 4), rep("Spleen_TEM", 4), rep("Liver_TEM", 2), rep("Liver_TRM", 2)))
)

Mackaydds <- DESeqDataSetFromMatrix(countData = Mackaycounts_clean[,Mtarget_columns],
                                    colData = Mcol_data,
                                    design = ~ condition)
Mackaydds <- DESeq(Mackaydds)

# Helper function to extract upregulated genes
get_up_genes <- function(dds_obj, contrast_vec) {
  res <- results(dds_obj, contrast = contrast_vec) %>% as.data.frame() %>% na.omit()
  rownames(res[which(res$log2FoldChange > 1 & res$padj < 0.05), ])
}

LiverTRM_genes  <- get_up_genes(Mackaydds, c("condition", "Liver_TRM", "Liver_TEM"))
LiverTEM_genes  <- get_up_genes(Mackaydds, c("condition", "Liver_TEM", "Spleen_TEM"))
SpleenTEM_genes <- get_up_genes(Mackaydds, c("condition", "Spleen_TEM", "Spleen_TCM"))
SpleenTCM_genes <- get_up_genes(Mackaydds, c("condition", "Spleen_TCM", "Spleen_TEM"))

writeLines(c(
  paste(c("MACKAY_LIVER_TRM_UP", "Upregulated_in_Liver_TRM_vs_TEM", LiverTRM_genes), collapse = "\t"),
  paste(c("MACKAY_LIVER_TEM_UP", "Upregulated_in_Liver_TEM_vs_Spleen_TEM", LiverTEM_genes), collapse = "\t"),
  paste(c("MACKAY_SPLEEN_TEM_UP", "Upregulated_in_Spleen_TEM_vs_Spleen_TCM", SpleenTEM_genes), collapse = "\t"),
  paste(c("MACKAY_SPLEEN_TCM_UP", "Upregulated_in_Spleen_TCM_vs_Spleen_TEM", SpleenTCM_genes), collapse = "\t")
), here("results", "tables", "Mackay_Four_Populations.gmt"))

################################################################################
# 3. GSEA & VISUALIZATIONS
################################################################################

# Extract specific ranked vectors directly from our pre-processed Module 2 list
# These replace res19, res22, res23, and res24 from the original script
get_ranked_vector <- function(contrast_name) {
  df <- mapped_res_list[[contrast_name]]
  # Ensure we have single mapping per symbol for ranking
  df_clean <- df %>% group_by(SYMBOL) %>% slice_max(abs(stat), n=1) %>% ungroup()
  setNames(df_clean$stat, df_clean$SYMBOL)
}

multi_sample_list <- list(
  "Ctrl_12hr_NoT vs Ctrl_12hr_T"   = get_ranked_vector("Ctrl_12hr_NoT_vs_Ctrl_12hr_T"),
  "Ctrl_12hr_NoT vs EKO_12hr_NoT"  = get_ranked_vector("Ctrl_12hr_NoT_vs_EKO_12hr_NoT"),
  "Ctrl_12hr_NoT vs S34_12hr_NoT"  = get_ranked_vector("Ctrl_12hr_NoT_vs_S34_12hr_NoT"),
  "Ctrl_12hr_NoT vs S4_12hr_NoT"   = get_ranked_vector("Ctrl_12hr_NoT_vs_S4_12hr_NoT")
)

combined_pathways <- list(
  "Wakim_Brain_TRM"   = brain_up, 
  "Mackay_Liver_TRM"  = LiverTRM_genes,
  "Mackay_Spleen_TCM" = SpleenTCM_genes,
  "Mackay_Spleen_TEM" = SpleenTEM_genes
)

# A. Multi-Sample Faceted Curves
multi_sample_curves <- map_df(names(combined_pathways), function(pw_name) {
  map_df(names(multi_sample_list), function(samp_name) {
    gsea_res <- tryCatch({
      plotEnrichmentData(combined_pathways[[pw_name]], multi_sample_list[[samp_name]])
    }, error = function(e) return(NULL))
    
    if (is.null(gsea_res)) return(NULL)
    
    tibble(
      rank    = as.numeric(gsea_res$curve$rank),
      ES      = as.numeric(gsea_res$curve$ES),
      Pathway = pw_name,
      Sample  = samp_name
    )
  })
})

# Save Faceted Plot
p_faceted <- ggplot(multi_sample_curves, aes(x = rank, y = ES, color = Sample)) +
  geom_line(linewidth = 1) +  # Updated from size to linewidth
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~Pathway, scales = "free_y") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "grey80", fill = NA),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  labs(title = "Multi-Sample GSEA Comparison", x = "Rank in Gene List", y = "Enrichment Score (ES)")

ggsave(here("results", "figures", "GSEA_Faceted_Comparison.png"), plot = p_faceted, width = 10, height = 8, dpi = 300)

# Save Separate Plots per Pathway
pathway_names <- unique(multi_sample_curves$Pathway)
for (pw in pathway_names) {
  p_single <- multi_sample_curves %>% 
    filter(Pathway == pw) %>%
    ggplot(aes(x = rank, y = ES, color = Sample)) +
    geom_line(linewidth = 1.2) + # Updated from size to linewidth
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    theme_minimal() +
    labs(title = paste("GSEA Profile:", pw), subtitle = "Overlaid Samples Comparison",
         x = "Rank in Gene List", y = "Enrichment Score (ES)") +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 14))
  
  ggsave(here("results", "figures", paste0("GSEA_Plot_", pw, ".png")), plot = p_single, width = 10, height = 6, dpi = 300)
}

# B. Smad4 vs Eomes Processing
shared_data <- read.csv(here("data", "raw_data", "shared_genes.csv"))

res_smad4_vector <- shared_data %>%
  mutate(stat = sign(lfc_S4) * -log10(pval_S4 + 1e-300)) %>%
  filter(!is.na(stat), !is.na(SYMBOL)) %>%
  { setNames(.$stat, .$SYMBOL) }

res_eomes_vector <- shared_data %>%
  mutate(stat = sign(lfc_Eomes) * -log10(pval_Eomes + 1e-300)) %>%
  filter(!is.na(stat), !is.na(SYMBOL)) %>%
  { setNames(.$stat, .$SYMBOL) }

smad4_eomes_list <- list(
  "Smad4_Controlled" = res_smad4_vector,
  "Eomes_Controlled" = res_eomes_vector
)

smad_eomes_curves <- map_df(names(combined_pathways), function(pw_name) {
  map_df(names(smad4_eomes_list), function(samp_name) {
    gsea_res <- tryCatch({
      plotEnrichmentData(combined_pathways[[pw_name]], smad4_eomes_list[[samp_name]])
    }, error = function(e) return(NULL))
    
    if (is.null(gsea_res) || is.null(gsea_res$curve)) return(NULL)
    
    tibble(
      rank    = as.numeric(gsea_res$curve$rank),
      ES      = as.numeric(gsea_res$curve$ES),
      Pathway = pw_name,
      Factor  = samp_name
    )
  })
})

for (pw in names(combined_pathways)) {
  plot_data_sub <- smad_eomes_curves %>% filter(Pathway == pw)
  if (nrow(plot_data_sub) == 0) next
  
  q <- ggplot(plot_data_sub, aes(x = rank, y = ES, color = Factor)) +
    geom_line(linewidth = 1.5) + # Updated from size to linewidth
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Smad4_Controlled" = "#D55E00", "Eomes_Controlled" = "#0072B2")) +
    theme_minimal() +
    labs(title = paste("Enrichment Profile:", pw), subtitle = "Overlay: Smad4 vs. Eomes Ranking",
         x = "Rank (Shared Genes)", y = "Enrichment Score (ES)") +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  
  ggsave(here("results", "figures", paste0("GSEA_Plot_shared_", pw, ".png")), plot = q, width = 10, height = 6, dpi = 300)
}

# C. Final Statistical Report
gsea_stats_report <- map_df(names(multi_sample_list), function(samp_name) {
  res <- fgsea(pathways = combined_pathways, 
               stats = multi_sample_list[[samp_name]],
               minSize = 10, maxSize = 500)
  res %>% mutate(Sample = samp_name) %>% as_tibble()
})

final_report <- gsea_stats_report %>%
  dplyr::select(Sample, Pathway = pathway, Size = size, ES = ES, NES = NES, 
                P_Value = pval, FDR_q_val = padj, leadingEdge) %>%
  arrange(Sample, desc(NES))

# Remove leadingEdge vector to write successfully to CSV
write.csv(final_report %>% dplyr::select(-leadingEdge), 
          here("results", "tables", "GSEA_MackayandWakim_Summary.csv"), 
          row.names = FALSE)