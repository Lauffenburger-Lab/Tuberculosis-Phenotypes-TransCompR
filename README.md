## Cross-species transcriptomics translation reveals a role for the unfolded protein response in Mycobacterium tuberculosis infection

Github repository for the study: 

Cross-species transcriptomics translation reveals a role for the unfolded protein response in Mycobacterium tuberculosis infection

Krista M. Pullen1*, Ryan Finethy2*, Seung-Hyun B. Ko1, Charlotte J. Reames2, Christopher M. Sassetti2^, Douglas A. Lauffenburger1^

(1) Department of Biological Engineering, Massachusetts Institute of Technology, Cambridge, MA

(2) Department of Microbiology and Physiological Systems, UMass Chan Medical School, Worcester, MA

(*) These authors contributed equally

(^) Co-corresponding authors

Douglas A Lauffenburger, 77 Massachusetts Avenue, Cambridge MA 02139, 617-252-1629 Email:  lauffen@mit.edu 

Christopher M. Sassetti, 368 Plantation Street, Worcester, MA 01605, 508-856-3678 
Email: christopher.sassetti@umassmed.edu

doi: 10.1038/s41540-024-00487-6

This repository is administered by @krista-pullen. For questions contact kpullen@alum.mit.edu

Numerous studies have identified similarities in blood transcriptomic signatures of tuberculosis (TB) phenotypes between mice and humans, including type 1 interferon production and innate immune cell activation. However, murine infection pathophysiology is distinct from human disease. We hypothesized that this is partly due to differences in the relative importance of biological pathways across species. To address this animal-to-human gap, we applied a systems modeling framework, Translatable Components Regression, to identify the axes of variation in the preclinical data most relevant to human TB disease state. Among the pathways our cross-species model pinpointed as highly predictive of human TB phenotype was the infection-induced unfolded protein response. To validate this mechanism, we confirmed that this cellular stress pathway modulates immune functions in Mycobacterium tuberculosis-infected mouse macrophages. Our work demonstrates how systems-level computational models enhance the value of animal studies for elucidating complex human pathophysiology.

# This repository contains: 

Code 

In the code folder, there is an R script for each figure that recreates the analyses and subplots from the manuscript. Figures without code do not include novel analyses (differential expression analyses for Figure S1) or utilize subscription-based software (Ingenuity Pathway Analysis and Prism for Figure 5-6 and S6). 

Data 

For the cross-species TransComp-R model, this study utilized publicly available mouse (GSE137092) and human (GSE107991 and GSE107992) bulk blood transcriptomics datasets. We downloaded the raw FASTQ files from NCBI GEO and performed STAR alignment (version 2.7.1a) to the Mm10 mouse and Hg38 human reference genomes, respectively. RSEM (version 1.3.1) was used to generate BAM files with gene and isoform expression level estimates from paired-end reads. The R package Tximport (version 1.16.1) was then used to generate count matrices from the gene expression level outputs from RSEM. Size factors were calculated and genes with fewer than 10 counts across the dataset were removed. Differential expression analysis was performed on the human data using DESeq2 (version 1.28.1) with Benjamini-Hochberg FDR correction. Differentially expressed genes (DEGs) were identified using a significance threshold of FDR<0.1 and abs(log2FoldChange)>0.5 to identify DEGs between latently and actively infected TB patients in either human cohort. DEGs were pooled across cohorts to create a single list of human DEGs. The DESeq2 variance-stabilizing transformation was applied to the count matrices for each geographic cohort separately, followed by z-score normalization. Finally, the normalized human counts matrices were aggregated and filtered to only include the DEGs. This aggregated human dataset was filtered for one-to-one mouse-human gene homologs identified by the Bioconductor tools biomaRt and homologene (versions 2.44.0 and 1.4.68.19.3.27). Genes where more than one gene paired to a single gene of the opposite species were removed. The normalized gene expression matrices including matched orthologs for each species can be found in the data folder of this repository (mouse_normalized_orthologs.csv and human_normalized_orthologs.csv), along with metadata files (mouse_metadata.csv and human_metadata.csv), and a dictionary of the matched orthologs (ortholog_dictionary.csv). The normalized but unfiltered gene expression matrices for each species are also in the data folder (mouse_normalized.csv and human_normalized.csv) but these datasets were only utilized in this study to generate Figure 1B. 
Note: the results of ortholog identification will vary slightly depending on the version of the biomaRt package and the ensembl database called (for the manuscript, we utilized ensembl archive release 104 – May 2021).

# R Dependencies: 

R version 4.3.2

biomaRt (version 2.44.0) 

caret (version 6.0-94)

corto (version 1.2.4)

dplyr (version 1.1.4)

effsize (version 0.8.1)

ggExtra (version 0.10.1)

ggplot2 (version 3.5.1)

ggpubr (version 0.6.0)

ggsignif (version 0.6.4)

pals (version 1.9)

pROC (version 1.18.5)

reshape2 (version 1.4.4)

TBSignatureProfiler (version 1.14.0)

