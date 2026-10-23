#!/bin/bash
#SBATCH --job-name=mint_cfg
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=Log/mint_cfg_%j.out
set -euo pipefail
cd <MINTCNA_ROOT>/ATAC_RNA
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
# Score every MintCNA configuration through the SAME projection the headline
# table uses (Code/eval/pearson_vs_truth.R), not a separate script, so the
# numbers are comparable to the benchmark rows rather than merely similar.
for TAG in a010e020 a010e025 a005e030; do
  echo "################ MINT_TAG=$TAG"
  MINT_TAG=$TAG CELLSET=own Rscript Code/eval/pearson_vs_truth.R P111T3 P110T1 P104T1 \
    2>/dev/null | grep -E "MintCNA|care_id" | head -8
  cp Result/qc/pearson_vs_truth_own.tsv Result/qc/mintcfg_${TAG}.tsv
done
echo "################ combined"
Rscript -e '
suppressPackageStartupMessages(library(data.table))
r <- rbindlist(lapply(c("a010e020","a010e025","a005e030"), function(t){
  d <- fread(sprintf("Result/qc/mintcfg_%s.tsv", t))[method=="MintCNA"]
  d[, tag := t][, .(tag, care_id, r = round(r,3), r_common = round(r_common,3), n_cell)] }))
w <- dcast(r, tag ~ care_id, value.var = "r")
w[, mean := round(rowMeans(.SD, na.rm=TRUE),3), .SDcols = setdiff(names(w),"tag")]
print(w); cat("\n"); print(r[order(care_id, -r)])'
