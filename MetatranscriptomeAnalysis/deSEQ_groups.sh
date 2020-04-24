#!/bin/bash
#SBATCH --mem=200000
#SBATCH --time=7-0:0:0
#SBATCH --mail-user="rstevick@my.uri.edu"
#SBATCH --mail-type=END,FAIL
#SBATCH --output="slurm-%x-%j.out"

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

# STEP 6: R ANALYSIS
# Note: For R to properly identify files to compare/contrast, they must include
# the appropriate prefix (either "control_$file" or experimental_$file")!


module load R-bundle-Bioconductor/3.3-foss-2016b-R-3.3.1

#checked Rscript $R_DIR/run_DESeq_stats_groups.R \
#  -I $STEP_5/RefSeq_results/org_results \
#  -O RefSeq_org_DESeq_results.tab \
#  -R $STEP_2/raw_counts.txt
#checked Rscript $R_DIR/run_DESeq_stats_groups.R \
#  -I $STEP_5/RefSeq_results/func_results \
#  -O RefSeq_func_DESeq_results.tab \
#  -R $STEP_2/raw_counts.txt

checked Rscript $R_DIR/Subsystems_DESeq_stats_comparisons.R \
  -I $STEP_5/Subsystems_results \
  -O Subsystems_level-4_DESeq_results.tab -L 4 \
  -R $STEP_2/raw_counts.txt



echo "Master bash script finished running at: "; date
exit
####################################################################
