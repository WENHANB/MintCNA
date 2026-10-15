#!/bin/sh
#SBATCH --job-name=MintCNA         # Job name
#SBATCH --output=log/MintCNA_%A_%a.out
#SBATCH --error=log/MintCNA_%A_%a.err
#SBATCH --account=feifeixiao	                  # your own sponser or account from this class
#SBATCH --qos=feifeixiao                     # your own sponser or account from this class
#SBATCH --mail-type=ALL                    # Mail events
#SBATCH --mail-user=<EMAIL>     # Where to send mail
#SBATCH --nodes=1                          # Use one node
#SBATCH --ntasks=1                      # Number of MPI ranks
#SBATCH --cpus-per-task=8               # Number of cores per MPI rank 
#SBATCH --mem=64gb                          # Total memory limit
#SBATCH --time=3-00:00:00                     # Time limit hrs:min:sec
#SBATCH --array=1-3


module load R/4.4 
BASE_LIST=("c3_N100_100kb_l1e4_n0.05" "c3_N100_100kb_l1e4_c0.05_n0.05" "c3_N100_100kb_l1e4_c0.2_n0.05")

BASE="${BASE_LIST[$SLURM_ARRAY_TASK_ID-1]}"

echo "Step1: RNA simulation"

Rscript <MINTCNA_ROOT>/Simulation/code/MintCNA/simulate_rna_from_cnasim.R "$BASE"

echo "Step2: DNA preprocess"

Rscript <MINTCNA_ROOT>/Simulation/code/MintCNA/DNA_process_n0.05.R "$BASE"

echo "Step3: RNA preprocess"

Rscript <MINTCNA_ROOT>/Simulation/code/MintCNA/RNA_process.R "$BASE"

echo "Step4: DNA and RNA combination"

Rscript <MINTCNA_ROOT>/Simulation/code/MintCNA/DNA_RNA_combine.R "$BASE"

echo "Step5: MintCNA"

Rscript <MINTCNA_ROOT>/Simulation/code/MintCNA/MintCNA.R "$BASE"

