#!/bin/bash
#SBATCH --job-name=rescore
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=Log/rescore_%j.out
set -euo pipefail
cd <MINTCNA_ROOT>/ATAC_RNA
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"

# All defaults come from the script itself (EPI_TAG w100kb_mf20000, RIDDLER_DIR
# 1MB, INFERCNV_TAG default, MINT_TAG a010e020). Nothing is overridden here --
# overriding in the launcher is exactly how the tuned settings got in last time.
echo "################ CELLSET=own"
CELLSET=own    Rscript Code/eval/pearson_vs_truth.R
echo "################ CELLSET=common"
CELLSET=common Rscript Code/eval/pearson_vs_truth.R

echo "################ figures"
for S in P111T3 P110T1 P104T1; do
  Rscript Code/figures/fig_truth_vs_methods.R "$S" || echo "  figure failed: $S"
done
echo "################ done"
