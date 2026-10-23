#!/bin/bash
#SBATCH --job-name=score_cfg
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=Log/score_cfg_%j.out
set -euo pipefail
cd <MINTCNA_ROOT>/ATAC_RNA
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
Rscript Code/eval/score_mint_configs.R
