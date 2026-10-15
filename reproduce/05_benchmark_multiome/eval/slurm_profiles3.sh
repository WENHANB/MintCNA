#!/bin/bash
#SBATCH --job-name=prof3
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=Log/prof3_%j.out
set -euo pipefail
cd <MINTCNA_ROOT>/ATAC_RNA
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
# Regenerate the projected tracks with the CHOSEN MintCNA configuration. The
# config sweep left them on a005e030, which is the worst of the three.
MINT_TAG=a010e020 CELLSET=own Rscript Code/eval/pearson_vs_truth.R P111T3 P110T1 P104T1
