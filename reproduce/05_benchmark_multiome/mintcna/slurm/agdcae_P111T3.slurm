#!/bin/bash
#SBATCH --job-name=agdcae_P111T3
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Code/mintcna/slurm/agdcae_P111T3_%j.log

set -euo pipefail
REPO=<MINTCNA_ROOT>/MintCNA_git
WORK=<MINTCNA_ROOT>/ATAC_RNA
OUT=$WORK/Result/MintCNA/P111T3_atac

module load conda
source activate cae_cbam
# run_denoise.py imports ag_dcae from python/, and mintcna.sh sets this too
export PYTHONPATH="$REPO/python:${PYTHONPATH:-}"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd "$WORK"
python "$REPO/scripts/run_denoise.py" \
  --input     "$OUT/P111T3_l2r_dna.csv" \
  --output    "$OUT/P111T3_ae_dna.csv" \
  --latent    "$OUT/P111T3_latent_z.csv" \
  --loss_plot "$OUT/P111T3_loss.pdf" \
  --config    "$WORK/Code/mintcna/config_P111T3.yaml"

echo "=== denoise done ==="
ls -la "$OUT"
