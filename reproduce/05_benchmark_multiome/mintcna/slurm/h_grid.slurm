#!/bin/bash
#SBATCH --job-name=h_grid
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=96gb
#SBATCH --array=1-6
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Code/mintcna/slurm/hgrid_%A_%a.log

# h sweep. The scan is the expensive stage and is cached per h by scan_file(h),
# so each array task does ONE scan then sweeps alpha x eps off it in seconds.
# eps 0.15 is included because 0.20 was the smallest value ever tried and it won
# at every sample -- the optimum may be outside the grid we searched.
set -euo pipefail
COMBOS=("P111T3 10" "P111T3 15" "P110T1 10" "P110T1 15" "P104T1 10" "P104T1 15")
read CARE H <<< "${COMBOS[$((SLURM_ARRAY_TASK_ID-1))]}"

REPO=<MINTCNA_ROOT>/MintCNA_git
WORK=<MINTCNA_ROOT>/ATAC_RNA
cd "$WORK"
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
CFG=$WORK/Code/mintcna/config_${CARE}_multi.yaml
OUT=$WORK/Result/MintCNA/${CARE}_multi
echo "=== $CARE  h=$H ==="

echo "--- scan (h=$H) ---"
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage scan --h "$H"

echo "--- call sweep ---"
for A in 0.05 0.10; do
  for E in 0.15 0.20 0.25; do
    TAG="h${H}a$(echo $A | tr -d '.')e$(echo $E | tr -d '0.')"
    TAG="h${H}a$(printf '%03d' $(echo "$A*100" | bc | cut -d. -f1))e$(printf '%03d' $(echo "$E*100" | bc | cut -d. -f1))"
    Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
      --h "$H" --alpha "$A" --st 3 --eps "$E" --tag "$TAG"
    cp "$OUT/${CARE}_mintcna_results.RData" "$OUT/${CARE}_results_${TAG}.RData"
    echo "  wrote $TAG"
  done
done
echo "=== $CARE h=$H done ==="
