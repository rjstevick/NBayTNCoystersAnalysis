#!/bin/bash
#SBATCH --cpus-per-task=20
#SBATCH --time=4-0:0:0
#SBATCH --mail-user="rstevick@my.uri.edu"
#SBATCH --mail-type=END,FAIL
#SBATCH --output="slurm-%x-%j.out"


# This script imports and QCs the data using DADA2
echo "QIIME2 bash script started running at: "; date

cd /data3/marine_diseases_lab/rebecca/tnc2017_16Sseq
module load QIIME2/2018.4 R/3.4.1-foss-2016b-X11-20160819
module list


# Put metadata file name here
METADATA="TNC_16S_Metadataall.txt"

# Import data into QIIME
# Paired-end, based on sample-manifest.csv
 qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path sample-manifest_TNC16Sall.csv --output-path TNC-paired-end-sequences.qza \
  --source-format PairedEndFastqManifestPhred33

# QC using dada2
qiime dada2 denoise-paired --verbose --i-demultiplexed-seqs TNC-paired-end-sequences.qza \
  --p-trunc-len-r 75 --p-trunc-len-f 75 \
  --p-trim-left-r 19 --p-trim-left-f 19 \
  --o-table table-notrim.qza \
  --o-representative-sequences rep-seqs-notrim.qza \
  --o-denoising-stats denoising-stats-notrim.qza \
  --p-n-reads-learn 100000000
  --p-n-threads 20

# Summarize feature table and sequences
qiime metadata tabulate \
  --m-input-file denoising-stats-notrim.qza \
  --o-visualization denoising-stats-notrim.qzv
qiime feature-table summarize \
  --i-table table-notrim.qza \
  --o-visualization table-notrim.qzv \
  --m-sample-metadata-file $METADATA
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-notrim.qza \
  --o-visualization rep-seqs-notrim.qzv

# This script assigns taxonomy based on the imported database (see X_qiime2gettaxonomydb.sh)

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs-notrim.qza \
  --o-classification taxonomy-notrim.qza
qiime metadata tabulate \
  --m-input-file taxonomy-notrim.qza \
  --o-visualization taxonomy-notrim.qzv
qiime taxa barplot \
  --i-table table-notrim.qza \
  --i-taxonomy taxonomy-notrim.qza \
  --m-metadata-file $METADATA \
  --o-visualization taxa-bar-plots-notrim.qzv

# This script calculated phylogenetic trees for the data

# align and mask sequences
qiime alignment mafft \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza
qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza

# calculate tree
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# calculate overall diversity
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 95 \
  --m-metadata-file $METADATA \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA \
  --m-metadata-column Station \
  --o-visualization core-metrics-results/unweighted-unifrac-station-significance.qzv \
  --p-pairwise
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA  \
  --m-metadata-column SampleType \
  --o-visualization core-metrics-results/unweighted-unifrac-group-significance.qzv \
  --p-pairwise

# This script calculates the rarefaction curve for the data
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 600 \
  --m-metadata-file $METADATA \
  --o-visualization alpha-rarefaction.qzv
  
  
  
  
  
echo "END $(date)"
