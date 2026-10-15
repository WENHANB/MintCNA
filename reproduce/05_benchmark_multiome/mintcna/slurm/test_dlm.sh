#!/bin/bash
#SBATCH --job-name=test_dlm
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --array=1-2
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Log/test_dlm_%A_%a.out
set -euo pipefail
cd <MINTCNA_ROOT>/ATAC_RNA
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
S=(P111T3 P110T1); Rscript Code/mintcna/test_dlm.R "${S[$((SLURM_ARRAY_TASK_ID-1))]}"
