#!/bin/bash
#SBATCH --job-name=refresh
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96gb
#SBATCH --array=1-20
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Code/mintcna/slurm/refresh_%A_%a.log

# The multivariate null simulated its Gaussian with the SDs laid out along the
# genome instead of across modalities (rep(sd_row, each = n) against a
# column-major matrix), which made tau 14-18% too strict. The scan SCORES are
# unaffected -- the null only supplies the threshold -- so this re-simulates the
# null in place and re-thresholds, rather than re-scanning.
#
# Both alphas are called: under the corrected null tau(0.05) = 0.595 is close to
# what tau(0.10) = 0.585 used to be, so alpha = 0.05 now reproduces the old
# operating point and alpha = 0.10 is genuinely looser than anything tried.
set -euo pipefail
S=(P59T2 P59T3 P100T1 P100T3 P103T1 P103T2 P104T1 P104T2 P105T1 P105T3
   P106T1 P106T2 P108T1 P108T2 P109T1 P109T2 P110T1 P110T2 P111T2 P111T3)
CARE=${S[$((SLURM_ARRAY_TASK_ID-1))]}
REPO=<MINTCNA_ROOT>/MintCNA_git
WORK=<MINTCNA_ROOT>/ATAC_RNA
cd "$WORK"
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
CFG=$WORK/Code/mintcna/config_${CARE}_multi.yaml
OUT=$WORK/Result/MintCNA/${CARE}_multi

Rscript Code/mintcna/refresh_null.R "$CARE" multi
for A in 0.05 0.10; do
  AA=$(printf '%03d' $(echo "$A*100" | bc | cut -d. -f1))
  TAG="fixnull_a${AA}e020"
  Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
    --alpha "$A" --st 3 --eps 0.20 --tag "$TAG" >/dev/null
  cp "$OUT/${CARE}_mintcna_results.RData" "$OUT/${CARE}_results_${TAG}.RData"
  Rscript -e "e<-new.env(); load('$OUT/${CARE}_results_${TAG}.RData', e); stopifnot(nrow(get('cna_state_matrix',e))>0)"
  echo "  wrote $TAG"
done
echo "=== $CARE done ==="
