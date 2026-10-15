#!/bin/bash
#SBATCH --job-name=epiSweep
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --array=1-6
#SBATCH --cpus-per-task=4
#SBATCH --mem=128gb
#SBATCH --time=72:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Logs/epiAneufinder/sweep_%A_%a.out

# Parameter search on ONE sample before spending 20x on any setting.
# P105T1 is the test case: its truth carries an unambiguous whole-arm 1p loss
# (log2 -0.72) that the published default misses in 99% of nuclei.
#
# Only two knobs can move this. The |z| < 1 collapse in assign_gainloss() is
# hard-coded, and uq/lq only trim within a segment, so neither is reachable.
#   windowSize : 1e5 (docs) vs 1e6 (roxygen, "Reccomended for sparse data")
#   k          : breakpoints per chromosome per cell; segments <= k+1. With 2
#                segments both z-scores are +/-0.707 and collapse to normal, so a
#                clean whole-arm event is structurally uncallable. More segments
#                may break the degeneracy.
set -u
export TMPDIR=$HOME/rtmp; mkdir -p "$TMPDIR"
module load R/4.4
ROOT=<MINTCNA_ROOT>/ATAC_RNA
CARE=P105T1
MF=5000

case $SLURM_ARRAY_TASK_ID in
  1) W=100000;  K=4  ;;   # published default window, default k  (baseline)
  2) W=100000;  K=8  ;;
  3) W=100000;  K=12 ;;
  4) W=1000000; K=4  ;;   # author's sparse-data window
  5) W=1000000; K=8  ;;
  6) W=1000000; K=12 ;;
esac
TAG="sweep_w$((W/1000))kb_k${K}"
echo "task $SLURM_ARRAY_TASK_ID -> $CARE window=$W k=$K tag=$TAG"
cd "$ROOT"
Rscript Code/benchmark/run_epiAneufinder.R "$CARE" "$W" "$MF" "$TAG" "$K"
