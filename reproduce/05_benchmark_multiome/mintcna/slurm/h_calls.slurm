#!/bin/bash
#SBATCH --job-name=h_calls
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --array=1-6
#SBATCH --output=<REPRODUCE_DIR>/05_benchmark_multiome/mintcna/slurm/hcalls_%A_%a.log

# Call-only sweep. Job 39710954 completed every SCAN (the expensive stage, cached
# as scan_h<h>.rds) but its call loops died when run_mintcna.R was edited under
# them. Nothing here re-scans.
set -euo pipefail
COMBOS=("P111T3 10" "P111T3 15" "P110T1 10" "P110T1 15" "P104T1 10" "P104T1 15")
read CARE H <<< "${COMBOS[$((SLURM_ARRAY_TASK_ID-1))]}"
REPO=<REPRODUCE_DIR>/..
WORK=<MINTCNA_ROOT>/ATAC_RNA
cd "$WORK"
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
CFG=$WORK/Code/mintcna/config_${CARE}_multi.yaml
OUT=$WORK/Result/MintCNA/${CARE}_multi
[ -s "$OUT/${CARE}_scan_h${H}.rds" ] || { echo "no scan checkpoint for h=$H"; exit 1; }
echo "=== $CARE h=$H ==="
for A in 0.05 0.10; do
  for E in 0.15 0.20 0.25; do
    TAG="h${H}a$(printf '%03d' $(echo "$A*100" | bc | cut -d. -f1))e$(printf '%03d' $(echo "$E*100" | bc | cut -d. -f1))"
    if [ -s "$OUT/${CARE}_results_${TAG}.RData" ]; then echo "  skip $TAG"; continue; fi
    Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
      --h "$H" --alpha "$A" --st 3 --eps "$E" --tag "$TAG" >/dev/null
    cp "$OUT/${CARE}_mintcna_results.RData" "$OUT/${CARE}_results_${TAG}.RData"
    echo "  wrote $TAG"
  done
done
echo "=== done ==="
