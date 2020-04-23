# Analysis bash scripts, processed data files, and R scripts for *Functional plasticity in oyster gut microbiomes along a eutrophication gradient in an urbanized estuary*
## Code written by Rebecca Stevick, URI-GSO (contact: rstevick@gmail.com)
### Collaborators: Anton F. Post and Marta Gómez-Chiarri

### This repository contains the scripts, pre-processed sequencing data, and the R script files to reproduce the figures in the manuscript. The raw sequences generated for this study can be found in the NCBI Short Read Archive under BioProject no. PRJNA598635.

### To cite this work: 
Stevick, R. J. (2019). Oyster-Associated Microbial Community Dynamics (Doctoral dissertation, University of Rhode Island).

-------------------------------------------------------------------------------------------------------------------------------

## Contents:
### MetatranscriptomeAnalysis
This folder contains scripts used to perform modified SAMSA2 analysis on the metatranscriptomic data.

### 16SampliconAnalysis
This folder contains scripts used to perform QIIME2 analysis on the 16S rRNA amplicon data.

### TNCpaper_ScriptsData
This folder contains all the processed data files and scripts to reproduce the figures and statistics in the manuscript. Each script is named for the figures it creates. All statistical analyses and plotting commands are included. 
#### Function
1. DeSeq2 results for level-2 SEED functional annnotation of metatranscriptomes (each station vs. mean)
2. DeSeq2 results for level-3 SEED functional annnotation of metatranscriptomes (each station vs. mean)
3. DeSeq2 results for level-4 SEED functional annnotation of metatranscriptomes (each station vs. mean)
4. DeSeq2 results for SEED functional annnotation of metatranscriptomes (North vs. South)
#### Metadata
1. Folder containing mapping files to creat Figure 1
2. Metadata spreadsheet with oyster and environmental measurements
#### Taxonomy
1. 16S rRNA taxonomy results
2. LefSe analysis results for the 16S rRNA taxa at the order level
3. Metatranscriptomic annotation taxonomy results with RefSeq
4. Metatranscriptomic annotation taxonomy results with Kraken (_RefSeq is used in the publication since we obtained higher annotation rates_)


