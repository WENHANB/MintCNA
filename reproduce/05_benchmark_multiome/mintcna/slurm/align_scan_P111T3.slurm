#!/bin/bash
#SBATCH --job-name=mint_scan_P111T3
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=200gb
#SBATCH --output=<REPRODUCE_DIR>/05_benchmark_multiome/mintcna/slurm/mint_scan_P111T3_%j.log

# align (cheap) + scan (the expensive, alpha-independent half).
# Calling is NOT run here: it takes seconds off the scan checkpoint, so it is a
# separate command you can repeat for every alpha / eps you want to try.
set -euo pipefail
REPO=<REPRODUCE_DIR>/..
WORK=<MINTCNA_ROOT>/ATAC_RNA
CFG=$WORK/Code/mintcna/config_P111T3_multi.yaml

module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"   # R/4.4 defaults to the 4.3 lib
export OMP_NUM_THREADS=1        # mclapply forks; BLAS threads on top would oversubscribe
export MKL_NUM_THREADS=1

cd "$WORK"
echo "=== align ==="
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage align

echo "=== scan (mc_cores from config) ==="
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage scan

echo "=== done ==="
ls -la "$WORK/Result/MintCNA/P111T3_multi/"
