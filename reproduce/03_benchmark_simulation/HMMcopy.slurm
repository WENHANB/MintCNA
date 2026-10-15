#!/bin/sh
#SBATCH --job-name=HMMcopy         # Job name
#SBATCH --output=log/HMMcopy_%A_%a.out
#SBATCH --error=log/HMMcopy_%A_%a.err
#SBATCH --account=feifeixiao	                  # your own sponser or account from this class
#SBATCH --qos=feifeixiao                     # your own sponser or account from this class
#SBATCH --mail-type=ALL                    # Mail events
#SBATCH --mail-user=<EMAIL>     # Where to send mail
#SBATCH --nodes=1                          # Use one node
#SBATCH --ntasks=1                      # Number of MPI ranks
#SBATCH --cpus-per-task=4               # Number of cores per MPI rank 
#SBATCH --mem=32gb                          # Total memory limit
#SBATCH --time=3-00:00:00                     # Time limit hrs:min:sec
#SBATCH --array=1-1

module load R/4.4 

BASE_LIST=("c3_N100_100kb_l1e4")
BASE="${BASE_LIST[$SLURM_ARRAY_TASK_ID-1]}"

START=$(date +%s)

Rscript <MINTCNA_ROOT>/Simulation/code/Benchmark/HMMcopy_new.R "$BASE"

END=$(date +%s)
ELAPSED=$(( (END - START) / 60 ))
echo "$BASE,${ELAPSED}min,$(sacct -j $SLURM_JOB_ID --format=MaxRSS --noheader | head -1)"