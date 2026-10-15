#!/bin/sh
#SBATCH --job-name=FLCNA         # Job name
#SBATCH --output=log/FLCNA_%A_%a.out
#SBATCH --error=log/FLCNA_%A_%a.err
#SBATCH --account=feifeixiao	                  # your own sponser or account from this class
#SBATCH --qos=feifeixiao-b                     # your own sponser or account from this class
#SBATCH --mail-type=ALL                    # Mail events
#SBATCH --mail-user=<EMAIL>     # Where to send mail
#SBATCH --nodes=1                          # Use one node
#SBATCH --ntasks=1                      # Number of MPI ranks
#SBATCH --cpus-per-task=8               # Number of cores per MPI rank 
#SBATCH --mem=64gb                          # Total memory limit
#SBATCH --time=3-00:00:00                     # Time limit hrs:min:sec
#SBATCH --array=1-2

module load R/4.4 

BASE_LIST=("c4_N250_100kb_l1e4" "c6_N500_100kb_l1e4")
BASE="${BASE_LIST[$SLURM_ARRAY_TASK_ID-1]}"

Rscript <REPRODUCE_DIR>/03_benchmark_simulation/FLCNA.R "$BASE"