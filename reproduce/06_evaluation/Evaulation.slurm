#!/bin/sh
#SBATCH --job-name=Evaulation         # Job name
#SBATCH --output=log/Evaulation_%A_%a.out
#SBATCH --error=log/Evaulation_%A_%a.err
#SBATCH --account=feifeixiao	                  # your own sponser or account from this class
#SBATCH --qos=feifeixiao                     # your own sponser or account from this class
#SBATCH --mail-type=ALL                    # Mail events
#SBATCH --mail-user=<EMAIL>     # Where to send mail
#SBATCH --nodes=1                          # Use one node
#SBATCH --ntasks=1                      # Number of MPI ranks
#SBATCH --cpus-per-task=2               # Number of cores per MPI rank 
#SBATCH --mem=16gb                          # Total memory limit
#SBATCH --time=3-00:00:00                     # Time limit hrs:min:sec
#SBATCH --array=1-9

module load R/4.4 

BASE_LIST=("c3_N100_100kb_l1e4_n0.05" "c3_N100_100kb_l1e4_c0.05_n0.05" "c3_N100_100kb_l1e4_c0.2_n0.05" 
"c3_N100_100kb_l1e4_n0.1" "c3_N100_100kb_l1e4_c0.05_n0.1" "c3_N100_100kb_l1e4_c0.2_n0.1"
"c3_N100_100kb_l1e4" "c3_N100_100kb_l1e4_c0.05" "c3_N100_100kb_l1e4_c0.2")
#BASE_LIST=("c3_N200_100kb_l1e4")
BASE="${BASE_LIST[$SLURM_ARRAY_TASK_ID-1]}"

Rscript <REPRODUCE_DIR>/06_evaluation/Evaulation_all.R "$BASE"