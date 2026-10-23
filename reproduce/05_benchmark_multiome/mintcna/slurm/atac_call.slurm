#!/bin/bash
#SBATCH --job-name=atac_call
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96gb
#SBATCH --array=1-20
#SBATCH --output=<REPRODUCE_DIR>/05_benchmark_multiome/mintcna/slurm/ataccall_%A_%a.log

# Call the cached ATAC-only scans.
#
# alpha = 0.05, the same nominal level as the multi arm. tau differs between the
# arms because each is calibrated against its own null, which is correct: alpha
# is the parameter, tau is what the data makes of it.
#
# eps = 0.28, NOT 0.20. eps is an absolute band on the fused segment mean z, and
# z is 1.18-1.51x wider without the RNA arm to average into it (measured on six
# samples; median ratio 1.39, matching eps 0.233-0.307). A single global 0.28 is
# used rather than a per-sample value so eps stays a fixed parameter rather than
# something estimated from the data; the per-sample table goes in the supplement.
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
Rscript <REPRODUCE_DIR>/05_benchmark_multiome/mintcna/make_config_atac.R "$CARE"          # picks up alpha = 0.05
CFG=$WORK/Code/mintcna/config_${CARE}_atac.yaml
OUT=$WORK/Result/MintCNA/${CARE}_atac
[ -s "$OUT/${CARE}_scan_h15.rds" ] || { echo "no scan for $CARE"; exit 1; }

# eps = 0.20 as well, so the supplement can show what the unscaled band would give
# 0.28 came from sd(z) scaling and from quantile matching on P111T3, and both
# overshot: it left ATAC-only calling 21.9% of the genome altered against 24.5%
# for ATAC+RNA. The neutral band covers the CENTRE of the distribution, whose
# shape differs from the tails that sd() weights, so a one-sample match does not
# transfer. 0.25 is the direct interpolation between the measured 21.9% (eps
# 0.28) and 33.1% (eps 0.20).
for spec in "0.25 h15a005e025"; do
  set -- $spec
  Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
    --alpha 0.05 --st 3 --eps "$1" --tag "$2" >/dev/null
  cp "$OUT/${CARE}_mintcna_results.RData" "$OUT/${CARE}_results_$2.RData"
  Rscript -e "e<-new.env(); load('$OUT/${CARE}_results_$2.RData', e); stopifnot(nrow(get('cna_state_matrix',e))>0)"
  echo "  wrote $2"
done
echo "=== $CARE done ==="
