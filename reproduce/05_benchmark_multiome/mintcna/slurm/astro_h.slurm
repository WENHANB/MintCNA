#!/bin/bash
#SBATCH --job-name=astro_h
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=96gb
#SBATCH --array=1-12
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Code/mintcna/slurm/astroh_%A_%a.log

# h was chosen on two Oligodendroglioma samples and one Astrocytoma, and MintCNA
# then lost to AtaCNV on all six Astrocytomas while winning 8 of 14
# Oligodendrogliomas. The hypothesis is that h = 15 (15 Mb) matches the whole-arm
# 1p/19q codeletion of oligodendroglioma and over-smooths the more focal
# astrocytoma events. Test it directly: smaller h on the astrocytomas.
#
# AG-DCAE output already exists for every sample, so this re-scans only.
set -euo pipefail
S=(P59T2 P59T3 P104T1 P104T2 P106T1 P106T2)
H=(5 10)
i=$((SLURM_ARRAY_TASK_ID-1)); CARE=${S[$((i/2))]}; HH=${H[$((i%2))]}

REPO=<MINTCNA_ROOT>/MintCNA_git
WORK=<MINTCNA_ROOT>/ATAC_RNA
cd "$WORK"
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
CFG=$WORK/Code/mintcna/config_${CARE}_multi.yaml
OUT=$WORK/Result/MintCNA/${CARE}_multi
echo "=== $CARE  h=$HH ==="

[ -s "$OUT/${CARE}_scan_h${HH}.rds" ] \
  && echo "--- scan h=$HH already cached ---" \
  || Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage scan --h "$HH"

for E in 0.20 0.25; do
  EE=$(printf '%03d' $(echo "$E*100" | bc | cut -d. -f1))
  TAG="h${HH}a010e${EE}"
  [ -s "$OUT/${CARE}_results_${TAG}.RData" ] && { echo "  skip $TAG"; continue; }
  Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
    --h "$HH" --alpha 0.10 --st 3 --eps "$E" --tag "$TAG" >/dev/null
  cp "$OUT/${CARE}_mintcna_results.RData" "$OUT/${CARE}_results_${TAG}.RData"
  Rscript -e "e<-new.env(); load('$OUT/${CARE}_results_${TAG}.RData', e); stopifnot(nrow(get('cna_state_matrix',e))>0)"
  echo "  wrote $TAG"
done
echo "=== done ==="
