#!/bin/bash
#SBATCH --job-name=astro_t
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=96gb
#SBATCH --array=1-6
#SBATCH --output=<REPRODUCE_DIR>/05_benchmark_multiome/mintcna/slurm/astrot_%A_%a.log

# MintCNA under-calls on the astrocytomas: recall 0.360 (gain) and 0.391 (loss)
# against 0.725 and 0.643 on the oligodendrogliomas, while PRECISION is higher
# there (0.72/0.79 vs 0.44/0.85). It is not calling the wrong things, it is
# calling too few. Three knobs move that, and they act at different stages:
#
#   h      scan window. A short event inside a 15 Mb window is diluted by the
#          neutral flanks and never produces a change-point at all.
#   alpha  tau is the (1-alpha) quantile of the null, so LARGER alpha means a
#          looser threshold and MORE breakpoints -- segments get finer.
#   eps    the TGMM neutral band. Smaller eps calls more of the segments that
#          were found non-neutral. It cannot recover an event that never became
#          its own segment.
#
# Sweeping all three separates them: if only h helps, the events are being lost
# in the scan; if eps helps at fixed h, they are being found and then dismissed.
# Scans for h in {5,10,15} already exist, so this is call-only.
set -euo pipefail
S=(P59T2 P59T3 P104T1 P104T2 P106T1 P106T2)
CARE=${S[$((SLURM_ARRAY_TASK_ID-1))]}
REPO=<REPRODUCE_DIR>/..
WORK=<MINTCNA_ROOT>/ATAC_RNA
cd "$WORK"
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
CFG=$WORK/Code/mintcna/config_${CARE}_multi.yaml
OUT=$WORK/Result/MintCNA/${CARE}_multi
echo "=== $CARE ==="
for HH in 5 10 15; do
  [ -s "$OUT/${CARE}_scan_h${HH}.rds" ] || { echo "  no scan h=$HH, skipping"; continue; }
  for A in 0.10 0.20; do
    for E in 0.10 0.15 0.20; do
      AA=$(printf '%03d' $(echo "$A*100" | bc | cut -d. -f1))
      EE=$(printf '%03d' $(echo "$E*100" | bc | cut -d. -f1))
      TAG="h${HH}a${AA}e${EE}"
      [ -s "$OUT/${CARE}_results_${TAG}.RData" ] && { echo "  skip $TAG"; continue; }
      Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
        --h "$HH" --alpha "$A" --st 3 --eps "$E" --tag "$TAG" >/dev/null 2>&1 \
        || { echo "  FAILED $TAG"; continue; }
      cp "$OUT/${CARE}_mintcna_results.RData" "$OUT/${CARE}_results_${TAG}.RData"
      echo "  wrote $TAG"
    done
  done
done
echo "=== done ==="
