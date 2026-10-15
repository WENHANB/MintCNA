#!/bin/bash
#SBATCH --job-name=dlmxh
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=72gb
#SBATCH --array=4-15
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Code/mintcna/slurm/cross_%A_%a.log

# h and smoothing are NOT separable -- both set an effective smoothing scale --
# so the two earlier grids found conditional optima, not the joint one. DLM was
# only ever tested at h=5, where it lost by 0.09.
#
# Two questions here:
#   1. does DLM still lose once h is large?   h in {10,15} x dV/dW in {100,1000}
#   2. has h turned over yet?                 h = 20, no DLM (15 was the edge)
set -euo pipefail
# h = 20 FIRST. The burst QoS runs about three tasks at a time, so array order is
# execution order, and the h = 20 answer is the one that gates the next decision
# (whether to lock h = 15 and run the remaining 17 samples).
COMBOS=(
  "P111T3 20 0"   "P110T1 20 0"   "P104T1 20 0"
  "P111T3 10 100" "P111T3 15 100" "P111T3 10 1000" "P111T3 15 1000"
  "P110T1 10 100" "P110T1 15 100" "P110T1 10 1000" "P110T1 15 1000"
  "P104T1 10 100" "P104T1 15 100" "P104T1 10 1000" "P104T1 15 1000")
read CARE H SM <<< "${COMBOS[$((SLURM_ARRAY_TASK_ID-1))]}"

REPO=<MINTCNA_ROOT>/MintCNA_git
WORK=<MINTCNA_ROOT>/ATAC_RNA
cd "$WORK"
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
CFG=$WORK/Code/mintcna/config_${CARE}_multi.yaml
OUT=$WORK/Result/MintCNA/${CARE}_multi
echo "=== $CARE  h=$H  smooth=$SM ==="

Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage scan \
  --h "$H" --smooth "$SM" --smooth-rna 0

# alpha 0.10 won at every h without smoothing; eps 0.20 and 0.25 tied at h=15.
for E in 0.20 0.25; do
  EE=$(printf '%03d' $(echo "$E*100" | bc | cut -d. -f1))
  if [ "$SM" -gt 0 ]; then TAG="h${H}sm${SM}a010e${EE}"; else TAG="h${H}a010e${EE}"; fi
  [ -s "$OUT/${CARE}_results_${TAG}.RData" ] && { echo "  skip $TAG"; continue; }
  Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
    --h "$H" --smooth "$SM" --smooth-rna 0 --alpha 0.10 --st 3 --eps "$E" --tag "$TAG" >/dev/null
  cp "$OUT/${CARE}_mintcna_results.RData" "$OUT/${CARE}_results_${TAG}.RData"
  # the parallel filesystem produced one 1.1 MB file of nul bytes in an earlier
  # run; verify the copy is loadable rather than discovering it at scoring time
  Rscript -e "e<-new.env(); load('$OUT/${CARE}_results_${TAG}.RData', e); stopifnot(nrow(get('cna_state_matrix',e))>0)"
  echo "  wrote $TAG"
done
echo "=== done ==="
