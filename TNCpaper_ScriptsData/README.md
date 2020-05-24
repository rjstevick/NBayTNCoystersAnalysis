## TNCpaper_ScriptsData
This folder contains all the processed data files from the analysis described in [MetatranscriptomeAnalysis](/MetatranscriptomeAnalysis) & [16SampliconAnalysis](/16SampliconAnalysis), and the scripts to reproduce the figures and statistics in the manuscript. Each script is named for the figures it creates.

**Software used:** R v3.4.1

#### Scripts
- **Figure1_2.R** - makes the map of sites and PCA of environmental data.
- **Figure3_S1.R** - 16S rRNA alpha- and beta-diversity and rarefaction plots.
- **Figure4_S1_S2_S4.R** - compares taxonomy generated from RefSeq/metatranscriptomes and 16S rRNA data with heatmaps and UpSetR plots (4 and S2). Metatranscriptome taxonomy rarefaction and beta-diversity plots (S1 and S4).
- **Figure5_6.R** - Metaranscriptomic functional data at pathway levels, subset for stress response (N vs S), nitrogen and phosphorus metabolism (each site vs mean) with barplots.
- **FigureS3_lefse.R** - 16S rRNA amplicon data at order level output from LefSe (Taxonomy/TNC_summaryLefSeResults.xlsx)
- **FigureS5-S8.R** - Metaranscriptomic functional data at gene level for stress response, nitrogen metabolism, phosphorus metabolism (heatmaps)

#### 1. Function/ (from [MetatranscriptomeAnalysis](/MetatranscriptomeAnalysis))
- DeSeq2 results for level-2 SEED functional annnotation of metatranscriptomes (each station vs. mean)
- DeSeq2 results for level-3 SEED functional annnotation of metatranscriptomes (each station vs. mean)
- DeSeq2 results for level-4 SEED functional annnotation of metatranscriptomes (each station vs. mean)
- DeSeq2 results for SEED functional annnotation of metatranscriptomes (North vs. South)

#### 2. Metadata/
- Folder containing mapping files to create Figure 1
- Metadata spreadsheet with oyster and environmental measurements

#### 3. Taxonomy/ (from [MetatranscriptomeAnalysis](/MetatranscriptomeAnalysis) & [16SampliconAnalysis](/16SampliconAnalysis))
- 16S rRNA taxonomy results
- LefSe analysis results for the 16S rRNA taxa at the order level
- Metatranscriptomic annotation taxonomy results with RefSeq

### To cite this work:
Stevick, R. J. (2019). Oyster-Associated Microbial Community Dynamics (Doctoral dissertation, University of Rhode Island).
