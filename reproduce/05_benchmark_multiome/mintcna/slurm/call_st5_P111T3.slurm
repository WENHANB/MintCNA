#!/bin/bash
#SBATCH --job-name=mint_st5_P111T3
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Code/mintcna/slurm/mint_st5_P111T3_%j.log
set -euo pipefail
REPO=<MINTCNA_ROOT>/MintCNA_git
WORK=<MINTCNA_ROOT>/ATAC_RNA
CFG=$WORK/Code/mintcna/config_P111T3_multi.yaml
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
cd "$WORK"
# st = 5 off the SAME cached scan, so the only thing that changes is the number
# of states the TGMM fits. eps sets the neutral band; the outer del/amp cut is
# fixed at 0.5 in gausianMixure.
for A in 0.05; do
  for E in 0.30 0.20; do
    TAG="st5a${A/./}e${E/./}"
    echo "=== st=5 alpha=$A eps=$E ($TAG) ==="
    Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
      --alpha "$A" --st 5 --eps "$E" --tag "$TAG"
    cp "$WORK/Result/MintCNA/P111T3_multi/P111T3_mintcna_results.RData" \
       "$WORK/Result/MintCNA/P111T3_multi/P111T3_results_${TAG}.RData"
  done
done
echo "=== done ==="
