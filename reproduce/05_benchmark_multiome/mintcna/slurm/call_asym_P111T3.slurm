#!/bin/bash
#SBATCH --job-name=mint_asym_P111T3
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24gb
#SBATCH --output=<REPRODUCE_DIR>/05_benchmark_multiome/mintcna/slurm/mint_asym_%j.log
set -euo pipefail
REPO=<REPRODUCE_DIR>/..
WORK=<MINTCNA_ROOT>/ATAC_RNA
CFG=$WORK/Code/mintcna/config_P111T3_multi.yaml
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
cd "$WORK"

# The decomposition showed the loss arm is fine (+0.592 of r, vs AtaCNV's +0.636)
# and the gain arm is where it is lost (+0.063 vs +0.131): MintCNA recovers only
# 41% of the true gain amplitude against AtaCNV's 102%. A symmetric neutral band
# is the suspect -- log2(3/2)=+0.585 against log2(1/2)=-1.000 means a symmetric
# band in log2 demands 1.7x more evidence for a gain. eps1 = 0.585*eps0 equalises
# them per copy; the sweep brackets that value.
#   eps0  eps1   note
#   0.30  0.300  symmetric baseline (= a010e030)
#   0.30  0.176  0.585 * eps0, the principled ratio
#   0.30  0.220
#   0.30  0.130
#   0.25  0.146  same ratio at a tighter loss band
#   0.20  0.117
run () {  # alpha eps0 eps1 tag
  echo "=== alpha=$1 eps0=$2 eps1=$3 ($4) ==="
  Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
    --alpha "$1" --st 3 --eps0 "$2" --eps1 "$3" --tag "$4"
  cp "$WORK/Result/MintCNA/P111T3_multi/P111T3_mintcna_results.RData" \
     "$WORK/Result/MintCNA/P111T3_multi/P111T3_results_$4.RData"
}
run 0.10 0.30 0.300 asym_e30_30
run 0.10 0.30 0.176 asym_e30_18
run 0.10 0.30 0.220 asym_e30_22
run 0.10 0.30 0.130 asym_e30_13
run 0.10 0.25 0.146 asym_e25_15
run 0.10 0.20 0.117 asym_e20_12
echo "=== done ==="
