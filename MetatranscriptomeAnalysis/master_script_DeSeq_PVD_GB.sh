#!/bin/bash
#SBATCH --mem=200000
#SBATCH --time=7-0:0:0
#SBATCH --mail-user="rstevick@my.uri.edu"
#SBATCH --mail-type=END,FAIL
#SBATCH --output="slurm-%x-%j.out"


# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if user is not submitting jobs to SLURM

####################################################################
#
# master_script.sh
# Created April 2017 by Sam Westreich, github.com/transcript
# This version modified February 21, 2018
#
####################################################################
#
# This script sets up and runs through ALL steps in the SAMSA pipeline
# before the analysis (which is done in R, likely in RStudio).  Each
# step is set up below.
#
# The steps are:
#   1.   Read cleaning with Trimmomatic
#   2.   Merging with PEAR, if applicable
#   3.   rRNA removal with SortMeRNA
#   4.   Annotation using DIAMOND (by default against the RefSeq database)
#   5.   Aggregation using analysis_counter.py
#   4.1  Annotation using DIAMOND against the Subsystems database
#   5.1  Aggregation using Subsystems-specific analysis counter.py
#   6.   Running R scripts to get DESeq statistical analysis.
#
# NOTE: BEFORE running this script, please run package_installation.bash
# and full_database_download.bash located at:
# https://github.com/transcript/samsa2/tree/master/setup in order to set
# up SAMSA2 dependencies and download full databases.
#
#######################################################################
#
echo -e "NOTE: Before running this script, please run package_installation.bash and full_database_download.bash located at https://github.com/transcript/samsa2/tree/master/setup in order to set up SAMSA2 dependencies.\n\n"
#
# VARIABLES - set starting location and starting files location pathways
#
source "/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA/bash_scripts/lib/common.sh"

export SAMSA=${SAMSA}:/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA

INPUT_DIR=/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA/input_files
OUT_DIR=/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA

STEP_1="$OUT_DIR/step_1_output"
STEP_2="$OUT_DIR/step_2_output"
STEP_3="$OUT_DIR/step_3_output"
STEP_4="$OUT_DIR/step_4_output"
STEP_5="$OUT_DIR/step_5_output"

if [[ -n "$USE_TINY" ]]; then
  # Diamond databases
  diamond_database="$SAMSA/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB"
  diamond_subsys_db="$SAMSA/setup_and_test/tiny_databases/subsys_db_TINY_24MB"
  # Aggregation databases
  RefSeq_db="$SAMSA/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB.fa"
  Subsys_db="$SAMSA/setup_and_test/tiny_databases/subsys_db_TINY_24MB.fa"
  # Use test output directories
  STEP_1="${STEP_1}_test"
  STEP_2="${STEP_2}_test"
  STEP_3="${STEP_3}_test"
  STEP_4="${STEP_4}_test"
  STEP_5="${STEP_5}_test"
else
  # Diamond databases
  diamond_database="/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA/full_databases/RefSeq_bac"
  diamond_subsys_db="/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA//full_databases/subsys_db"
  # Aggregation databases
  RefSeq_db="/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA/full_databases/RefSeq_bac.fa"
  Subsys_db="/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA/full_databases/subsys_db.fa"
fi

####################################################################
#
# STEP 1: CLEANING FILES WITH TRIMMOMATIC


####################################################################
#
# STEP 2: MERGING OF PAIRED-END FILES USING PEAR

####################################################################
#
# STEP 2.9: GETTING RAW SEQUENCES COUNTS
# Note: These are used later for statistical analysis.

####################################################################
#
# STEP 3: REMOVING RIBOSOMAL READS WITH SORTMERNA
# Note: this step assumes that the SortMeRNA databases are indexed.  If not,
# do that first (see the SortMeRNA user manual for details).


####################################################################
#
# STEP 4: ANNOTATING WITH DIAMOND AGAINST REFSEQ
# Note: this step assumes that the DIAMOND database is already built.  If not,
# do that first before running this step.


####################################################################
#
# STEP 5: AGGREGATING WITH ANALYSIS_COUNTER

####################################################################
#
# STEP 4.1: ANNOTATING WITH DIAMOND AGAINST SUBSYSTEMS


##################################################################
#
# STEP 5.1: PYTHON SUBSYSTEMS ANALYSIS COUNTER

##################################################################
#
# At this point, all the results files are ready for analysis using R.
# This next step performs basic DESeq2 analysis of the RefSeq organism, function,
# and Subsystems annotations.
#
# More complex R analyses may be performed using specific .sh analysis scripts.
#
# STEP 6: R ANALYSIS
# Note: For R to properly identify files to compare/contrast, they must include
# the appropriate prefix (either "control_$file" or experimental_$file")!


module load R-bundle-Bioconductor/3.3-foss-2016b-R-3.3.1


checked Rscript $R_DIR/Subsystems_pie_charts.R \
  -I $STEP_5/Subsystems_results \
  -O Subsystems_level-1_pie_results.tab -L 1 \
  -R $STEP_2/raw_counts.txt

checked Rscript $R_DIR/run_DESeq_stats.R \
  -I $STEP_5/RefSeq_results/org_results \
  -O RefSeq_org_DESeq_results.tab \
  -R $STEP_2/raw_counts.txt
checked Rscript $R_DIR/run_DESeq_stats.R \
  -I $STEP_5/RefSeq_results/func_results \
  -O RefSeq_func_DESeq_results.tab \
  -R $STEP_2/raw_counts.txt
checked Rscript $R_DIR/Subsystems_DESeq_stats.R \
  -I $STEP_5/Subsystems_results \
  -O Subsystems_level-1_DESeq_results.tab -L 1 \
  -R $STEP_2/raw_counts.txt


echo "Master bash script finished running at: "; date
exit
####################################################################
