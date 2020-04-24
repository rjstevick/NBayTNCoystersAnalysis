#!/bin/bash
#SBATCH --ntasks-per-node=5 #processor cores per node
#SBATCH --time=8-0:0:0
#SBATCH --mail-user="rstevick@my.uri.edu"
#SBATCH --mail-type=END,FAIL
#SBATCH --output="slurm-%x-%j.out"
#SBATCH --exclusive=user


INPUT_DIR=/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA/input_files
OUT_DIR=/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA

  # diamond databases
diamond_database="/net/fs03/data3/marine_diseases_lab/shared/NCBI_NRdb/ncbinrdiamond"
diamond_subsys_db="/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA/full_databases/subsys_db"
  # Aggregation databases
nr_db="/net/fs03/data3/marine_diseases_lab/shared/NCBI_NRdb/nr.gz"
Subsys_db="/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA/full_databases/subsys_db.fa"


module load DIAMOND/0.9.23-foss-2016b


echo "Now starting on DIAMOND org annotations at: "; date

for file in /net/fs03/data3/marine_diseases_lab/rebecca/tnc2017_sequencing/08_SAMSA/step_3_output/*ribodepleted.fastq
do
    shortname=`echo $file | awk -F "ribodepleted" '{print $1 "RefSeq_annotated"}'`
    echo "Now starting on " $file
    echo "Converting to " $shortname
    diamond blastx --db $diamond_database -q $file -a $file.nrannot -t ./ -k 1 --threads 5
    diamond view --daa $file.nrannot.daa -o $shortname -f tab --threads 5
done

echo "RefSeq DIAMOND annotations completed at: "; date

####################################################################

