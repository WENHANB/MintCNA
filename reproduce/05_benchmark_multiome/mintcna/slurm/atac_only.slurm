#!/bin/bash
#SBATCH --job-name=atac_only
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=96gb
#SBATCH --array=1-20
#SBATCH --output=<REPRODUCE_DIR>/05_benchmark_multiome/mintcna/slurm/ataconly_%A_%a.log

# The ablation: MintCNA on the ATAC arm alone, everything else identical.
# No preprocess and no AG-DCAE -- both are shared with the multi run, so the only
# difference between the two arms is whether RNA enters the fusion.
# DNA-only mode uses the univariate CUSUM (p = 1), which is faster than the
# multivariate scan, so 4 h is generous.
set -euo pipefail
S=(P59T2 P59T3 P100T1 P100T3 P103T1 P103T2 P104T1 P104T2 P105T1 P105T3
   P106T1 P106T2 P108T1 P108T2 P109T1 P109T2 P110T1 P110T2 P111T2 P111T3)
CARE=${S[$((SLURM_ARRAY_TASK_ID-1))]}
REPO=<REPRODUCE_DIR>/..
WORK=<MINTCNA_ROOT>/ATAC_RNA
cd "$WORK"
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "=== $CARE  ATAC-only  ($(date)) ==="

Rscript <REPRODUCE_DIR>/05_benchmark_multiome/mintcna/make_config_atac.R "$CARE"
CFG=$WORK/Code/mintcna/config_${CARE}_atac.yaml
OUT=$WORK/Result/MintCNA/${CARE}_atac
TAG=h15a010e020

# SCAN ONLY. alpha and eps are call-stage parameters and cost seconds, and both
# still have to be settled for this arm: the corrected null moved tau, and the
# fused segment mean z is 1.40x wider without the RNA arm (sd 0.394 against
# 0.281 on P111T3), so eps = 0.20 here would not mean what it means in multi
# mode. Scanning now and calling later keeps the expensive half off the critical
# path of that decision.
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage align
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage scan
[ -s "$OUT/${CARE}_scan_h15.rds" ] || { echo "scan checkpoint missing"; exit 1; }
echo "  scan cached: $OUT/${CARE}_scan_h15.rds"
echo "=== $CARE done ($(date)) ==="
