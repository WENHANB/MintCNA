#!/bin/bash
#SBATCH --job-name=mint_call_P111T3
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Code/mintcna/slurm/mint_call_P111T3_%j.log

# Threshold + TGMM off the cached scan. Sized from the measured run: 217 s wall,
# 1.3 GB peak, single-threaded -- so 2 cpu / 8 GB, not the 32/200 the scan needed.
set -euo pipefail
REPO=<MINTCNA_ROOT>/MintCNA_git
WORK=<MINTCNA_ROOT>/ATAC_RNA
CFG=$WORK/Code/mintcna/config_P111T3_multi.yaml
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
cd "$WORK"

for A in 0.05 0.10; do
  for E in 0.30 0.25 0.20; do
    TAG="a${A/./}e${E/./}"
    echo "=== alpha=$A eps=$E  ($TAG) ==="
    Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
      --alpha "$A" --st 3 --eps "$E" --tag "$TAG"
    cp "$WORK/Result/MintCNA/P111T3_multi/P111T3_mintcna_results.RData" \
       "$WORK/Result/MintCNA/P111T3_multi/P111T3_results_${TAG}.RData"
  done
done
echo "=== done ==="
