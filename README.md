# Analysis bash scripts, processed data files, and R scripts for *Functional plasticity in oyster gut microbiomes along a eutrophication gradient in an urbanized estuary*

## Code written by Rebecca Stevick, URI-GSO (contact: rjstevick (at) gmail.com)
### Collaborators: Anton F. Post and Marta Gómez-Chiarri

### *Study of wild oysters and their gut microbiomes along the estuarine gradient of Narragansett Bay, Rhode Island, USA with 16S rRNA amplicon and metatranscriptome data.*

This repository contains the scripts, pre-processed sequencing data, and the R script files to reproduce the figures in the manuscript. The raw sequences generated for this study can be found in the NCBI Short Read Archive under BioProject no. PRJNA598635 [(reviewer link here)](https://dataview.ncbi.nlm.nih.gov/object/PRJNA598635?reviewer=5g8ftj52399ut5ujmj7v7t0hh8). The corresponding accession numbers for each sample are detailed in [`PRJNA598635_NCBI_info.xlsx`](/PRJNA598635_NCBI_info.xlsx). This file also includes file names and environmental data for each 16S rRNA amplicon or metatranscriptomic sample.

### To cite this work:
Stevick, R. J. (2019). Oyster-Associated Microbial Community Dynamics (Doctoral dissertation, University of Rhode Island).

--------------------------------------------------------------------------------

# Contents
## [MetatranscriptomeAnalysis](/MetatranscriptomeAnalysis)
This folder contains scripts used to perform modified SAMSA2 analysis on the metatranscriptomic data.
#### 1. Bash Scripts
- Main modified SAMSA2 script (masterscript*.sh)
- Scripts to run Rscripts on the command line (deSEQ*.sh)
- Script to compute DIAMOND annotations against refSeq dbs (newannot.sh)
- Raw sequence counts and file inputs (rawcounts.txt and samples_transcriptomes.txt)
#### 2. R_scripts
- DeSeq2 files to run states for the susbsytems based on deSEQ_groups.sh or deSEQ_subsys.sh


## [16SampliconAnalysis](/16SampliconAnalysis)
This folder contains scripts used to perform QIIME2 analysis on the 16S rRNA amplicon data.
#### 1. Scripts
- Main QIIME2 script with all commands (00_qiime2_alltogether.sh)
- Script to export QIIME2 files (qiime2_exportcollapsedtable.sh)
#### 2. QIIME2 Results
- *core-metrics-results*: diversity metrics calculated by QIIME2
- *exported*: files exported from QIIME2, then imported to excel
- *rarefaction*: curves generated by QIIME2
- *taxonomy*: annotations and barplots generated by QIIME2
- all other qza and qzv files generated during the QIIME2 pipeline


## [TNCpaper_ScriptsData](/TNCpaper_ScriptsData)
This folder contains all the processed data files from the analysis described in [MetatranscriptomeAnalysis](/MetatranscriptomeAnalysis) & [16SampliconAnalysis](/16SampliconAnalysis), and the scripts to reproduce the figures and statistics in the manuscript. Each script is named for the figures it creates.
#### 1. Function (from [MetatranscriptomeAnalysis](/MetatranscriptomeAnalysis))
- DeSeq2 results for level-2 SEED functional annnotation of metatranscriptomes (each station vs. mean)
- DeSeq2 results for level-3 SEED functional annnotation of metatranscriptomes (each station vs. mean)
- DeSeq2 results for level-4 SEED functional annnotation of metatranscriptomes (each station vs. mean)
- DeSeq2 results for SEED functional annnotation of metatranscriptomes (North vs. South)
#### 2. Metadata
- Folder containing mapping files to create Figure 1
- Metadata spreadsheet with oyster and environmental measurements
#### 3. Taxonomy (from [MetatranscriptomeAnalysis](/MetatranscriptomeAnalysis) & [16SampliconAnalysis](/16SampliconAnalysis))
- 16S rRNA taxonomy results
- LefSe analysis results for the 16S rRNA taxa at the order level
- Metatranscriptomic annotation taxonomy results with RefSeq
#### Scripts
- **Figure1_2_Table1.R** - makes the map of sites and PCA of environmental data, calculates environmental correlation coeffs for Table 1.
- **Figure3_S1_S2_S6.R** - 16S rRNA percent abundance barplots (including controls) (S1), rarefaction plots (S2), and alpha- and beta-diversity (3). 16S core microbiome (S6). Table S2 of beta-diversity stats.
- **Figure4_S2_S3_S5_S6.R** - compares taxonomy generated from RefSeq/metatranscriptomes and 16S rRNA data with heatmaps and UpSetR plots (4 and S3). Metatranscriptome taxonomy rarefaction (S2) and beta-diversity plots (S5). Metatranscriptome core microbiome (S6).
- **Figure5_6.R** - Metaranscriptomic functional data at pathway levels, subset for stress response (N vs S), nitrogen and phosphorus metabolism (each site vs mean) with barplots.
- **FigureS4_lefse.R** - 16S rRNA amplicon data at order level output from LefSe (Taxonomy/TNC_summaryLefSeResults.xlsx)
- **FigureS7-S10.R** - Metaranscriptomic functional data at gene level for stress response, nitrogen metabolism, phosphorus metabolism (heatmaps)
