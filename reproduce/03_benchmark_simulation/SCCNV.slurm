#!/bin/sh
#SBATCH --job-name=SCCNV         # Job name
#SBATCH --output=log/SCCNV_%A_%a.out
#SBATCH --error=log/SCCNV_%A_%a.err
#SBATCH --account=feifeixiao	                  # your own sponser or account from this class
#SBATCH --qos=feifeixiao-b                     # your own sponser or account from this class
#SBATCH --mail-type=ALL                    # Mail events
#SBATCH --mail-user=<EMAIL>     # Where to send mail
#SBATCH --nodes=1                          # Use one node
#SBATCH --ntasks=1                      # Number of MPI ranks
#SBATCH --cpus-per-task=8               # Number of cores per MPI rank 
#SBATCH --mem=64gb                          # Total memory limit
#SBATCH --time=3-00:00:00                     # Time limit hrs:min:sec
#SBATCH --array=1-9

module load R/4.4 

BASE_LIST=(
  "c3_N100_100kb_l1e4"
  "c3_N100_100kb_l1e4_c0.2"
  "c3_N100_100kb_l1e4_c0.2_n0.1"
  "c3_N100_100kb_l1e4_c0.2_n0.05"
  "c3_N100_100kb_l1e4_c0.05"
  "c3_N100_100kb_l1e4_c0.05_n0.1"
  "c3_N100_100kb_l1e4_c0.05_n0.05"
  "c3_N100_100kb_l1e4_n0.1"
  "c3_N100_100kb_l1e4_n0.05"
)
BASE="${BASE_LIST[$SLURM_ARRAY_TASK_ID-1]}"

START=$(date +%s)

Rscript <REPRODUCE_DIR>/03_benchmark_simulation/SCCNV.R "$BASE"

END=$(date +%s)
ELAPSED=$(( (END - START) / 60 ))
echo "$BASE,${ELAPSED}min,$(sacct -j $SLURM_JOB_ID --format=MaxRSS --noheader | head -1)"