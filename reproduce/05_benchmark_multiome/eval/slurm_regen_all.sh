#!/bin/bash
#SBATCH --job-name=regen
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Log/regen_%j.out
set -euo pipefail
cd <MINTCNA_ROOT>/ATAC_RNA
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
# Everything downstream reads Result/qc/profiles/*.tsv.gz and
# pearson_vs_truth_own.tsv, both still holding the pre-null-fix MintCNA
# (h15a010e020). One pass regenerates them at CHOSEN_MINT = h15a005e020 and adds
# the MintCNA(ATAC) arm. The six benchmark methods take no input from MintCNA, so
# their `r` must come back unchanged -- that is the check that nothing else moved.
CELLSET=own    Rscript Code/eval/pearson_vs_truth.R
CELLSET=common Rscript Code/eval/pearson_vs_truth.R
