#!/bin/bash
#SBATCH --job-name=sm_grid
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=96gb
#SBATCH --array=1-9
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Code/mintcna/slurm/smgrid_%A_%a.log

# End-to-end test of the genome-direction smoother. The pre-segmentation matrix
# improves a lot (P111T3 pseudobulk r 0.358 -> 0.750) but the amplitude at real
# events drops 0.505 -> 0.374, and our change-point scan finds breakpoints BY
# amplitude. Only the final score settles whether the trade is worth it, so this
# runs scan + call, not the diagnostic.
#
# dV/dW is the only knob: the smoother sees just the ratio. 100 was the best of
# {10, 33, 100} and the trend had not turned, so the grid extends upward.
set -euo pipefail
COMBOS=("P111T3 100" "P111T3 300" "P111T3 1000"
        "P110T1 100" "P110T1 300" "P110T1 1000"
        "P104T1 100" "P104T1 300" "P104T1 1000")
read CARE SM <<< "${COMBOS[$((SLURM_ARRAY_TASK_ID-1))]}"

REPO=<MINTCNA_ROOT>/MintCNA_git
WORK=<MINTCNA_ROOT>/ATAC_RNA
cd "$WORK"
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
CFG=$WORK/Code/mintcna/config_${CARE}_multi.yaml
OUT=$WORK/Result/MintCNA/${CARE}_multi
echo "=== $CARE  smooth=$SM ==="

Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage scan --smooth "$SM" --smooth-rna 0

# alpha/eps held at the values chosen without smoothing, so the comparison is
# against the current best and not against a re-tuned variant.
for spec in "0.10 0.20 a010e020" "0.10 0.15 a010e015" "0.05 0.20 a005e020"; do
  set -- $spec
  TAG="sm${SM}$3"
  Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
    --smooth "$SM" --smooth-rna 0 --alpha "$1" --st 3 --eps "$2" --tag "$TAG"
  cp "$OUT/${CARE}_mintcna_results.RData" "$OUT/${CARE}_results_${TAG}.RData"
  echo "  wrote $TAG"
done
echo "=== done ==="
