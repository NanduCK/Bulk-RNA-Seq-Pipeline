# Bulk RNA-Seq Analysis Pipeline: Investigating the Smad4/Eomes Axis

This repository contains a reproducible, modular R pipeline for processing and analysing bulk RNA-sequencing data. Developed to investigate the transcriptional regulation of the Smad4/Eomes axis, this workflow handles everything from raw count ingestion to differential expression and custom Gene Set Enrichment Analysis (GSEA).

To benchmark our findings against established immunological signatures, this pipeline also automatically ingests, normalises, and compares our data against published external datasets (Wakim et al. and Mackay et al.) to identify Tissue-Resident Memory (TRM) and Effector Memory (TEM) profiles.

Additionally, the Bulk_RNASeq_Smad4_Eomes_axis.R file contains the code initially written for the project. I later split this into the three other scripts for the pipeline to ensure proper modularity and reproducibility. 

## 🗂️ Repository Structure

The project is organised to ensure reproducibility across different machines using the `here` package. 

```text
├── README.md                           # Project overview and instructions
├── scripts/                            # R pipeline modules
   ├── 01_process_counts.R             # Data ingestion and metadata generation
   ├── 02_run_deseq2.R                 # Differential expression and ID mapping
   ├── 03_gsea_analysis.R              # External data integration and visualisation
   └──Bulk_RNASeq_Smad4_Eomes_axis.R   # Initial Pipeline 

```
