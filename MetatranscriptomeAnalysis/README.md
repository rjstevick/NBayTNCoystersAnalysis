## MetatranscriptomeAnalysis
This folder contains scripts used to perform modified SAMSA2 analysis on the metatranscriptomic data. The raw sequences generated for this study can be found in the NCBI SRA under BioProject no. PRJNA598635. See [`PRJNA598635_NCBI_info.xlsx`](/PRJNA598635_NCBI_info.xlsx) for accession numbers for each sample.

#### 1. Bash Scripts
- Main modified SAMSA2 script (masterscript*.sh)
- Scripts to run Rscripts on the command line (deSEQ*.sh)
- Script to compute DIAMOND annotations against refSeq dbs (newannot.sh)
- Raw sequence counts and file inputs (rawcounts.txt and samples_transcriptomes.txt)

#### 2. R_scripts
- DeSeq2 files to run states for the susbsytems based on deSEQ_groups.sh or deSEQ_subsys.sh

### To cite this work:
Stevick, R.J., Post, A.F. & Gómez-Chiarri, M (2021). Functional plasticity in oyster gut microbiomes along a eutrophication gradient in an urbanized estuary. *Animal Microbiome* 3, 5. https://doi.org/10.1186/s42523-020-00066-0
