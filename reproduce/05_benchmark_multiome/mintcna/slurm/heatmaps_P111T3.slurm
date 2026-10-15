#!/bin/bash
#SBATCH --job-name=mint_hm_P111T3
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16gb
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Code/mintcna/slurm/mint_hm_%j.log
set -euo pipefail
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
cd <MINTCNA_ROOT>/ATAC_RNA
Rscript Code/figures/mintcna_heatmaps.R P111T3
