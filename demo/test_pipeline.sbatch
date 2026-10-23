#!/bin/bash
#SBATCH --job-name=mintcna_demo
#SBATCH --account=feifeixiao
#SBATCH --partition=hpg-default
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --output=/orange/feifeixiao/wenhan.bao/Multi-omics/MintCNA_git/demo/test_%j.out
#SBATCH --error=/orange/feifeixiao/wenhan.bao/Multi-omics/MintCNA_git/demo/test_%j.err

set -euo pipefail

REPO=/orange/feifeixiao/wenhan.bao/Multi-omics/MintCNA_git
cd "$REPO"

# R and pytorch modules conflict (both use apptainer), so each stage runs in
# its own subshell with only the module it needs.
export R_LIBS=/home/wenhan.bao/R/x86_64-pc-linux-gnu-library/4.4

echo "=========================================="
echo "MintCNA Demo Pipeline Test"
echo "Start: $(date)"
echo "Host:  $(hostname)"
echo "=========================================="

echo ""
echo ">>> Stage 1: R preprocessing"
(
  module load R/4.4
  Rscript scripts/run_mintcna.R --config demo/demo_config.yaml --stage preprocess
)

echo ""
echo ">>> Stage 2: Python AG-DCAE denoising"
(
  module load pytorch/2.8.0
  export PYTHONPATH="$REPO/MintCNA:${PYTHONPATH:-}"
  python scripts/run_denoise.py \
    --input demo/output/l2r_dna.csv \
    --output demo/output/ae_dna.csv \
    --latent demo/output/latent_z.csv \
    --loss_plot demo/output/loss_curve.pdf \
    --config demo/demo_config.yaml
)

echo ""
echo ">>> Stage 3: R alignment + segmentation + calling"
(
  module load R/4.4
  Rscript scripts/run_mintcna.R --config demo/demo_config.yaml --stage all_post_denoise
)

echo ""
echo "=========================================="
echo "End: $(date)"
echo "=========================================="
echo ""
echo "Input files (demo/input/):"
ls -la demo/input/
echo ""
echo "Output files (demo/output/):"
ls -la demo/output/
