#!/bin/bash
#SBATCH --job-name=mint_resume
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=96gb
#SBATCH --array=1-2
#SBATCH --output=<REPRODUCE_DIR>/05_benchmark_multiome/mintcna/slurm/mint_resume_%A_%a.log

# Steps 5-6 of full_pipeline.slurm only. Job 39636410 finished AG-DCAE for both
# samples (3h09 and 1h30) and then died on `Rscript: command not found`, so the
# denoised matrices and latents are already on disk and the training must not be
# repeated. No conda in this script at all.
set -euo pipefail
SAMPLES=(P110T1 P104T1)
CARE=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
REPO=<REPRODUCE_DIR>/..
WORK=<MINTCNA_ROOT>/ATAC_RNA
cd "$WORK"
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
CFG=$WORK/Code/mintcna/config_${CARE}_multi.yaml
OUT=$WORK/Result/MintCNA/${CARE}_multi
echo "=== $CARE ==="
for f in "$OUT/${CARE}_ae_dna.csv" "$OUT/${CARE}_l2r_rna.csv" "$CFG"; do
  [ -s "$f" ] || { echo "missing prerequisite: $f"; exit 1; }
done

echo "--- 5. align + scan ---"
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage align
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage scan

echo "--- 6. call ---"
for spec in "0.10 0.20 a010e020" "0.10 0.25 a010e025" "0.05 0.30 a005e030"; do
  set -- $spec
  Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
    --alpha "$1" --st 3 --eps "$2" --tag "$3"
  cp "$OUT/${CARE}_mintcna_results.RData" "$OUT/${CARE}_results_$3.RData"
done
echo "=== $CARE done ==="
