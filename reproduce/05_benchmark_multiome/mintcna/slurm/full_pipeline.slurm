#!/bin/bash
#SBATCH --job-name=mint_full
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=96gb
#SBATCH --array=1-2
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Code/mintcna/slurm/mint_full_%A_%a.log

# Whole MintCNA chain for one sample: bin -> input -> preprocess -> AG-DCAE
# -> align -> scan -> call. Sized from the P111T3 measurements: AG-DCAE peaked at
# 2.9 GB, the scan at 21 GB, so 96 GB has plenty of headroom; 32 cores because
# the scan is the part that parallelises.
set -euo pipefail
SAMPLES=(P110T1 P104T1)
CARE=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
echo "=== $CARE ==="

REPO=<MINTCNA_ROOT>/MintCNA_git
WORK=<MINTCNA_ROOT>/ATAC_RNA
cd "$WORK"
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1     # mclapply forks; BLAS on top would oversubscribe

echo "--- 1. bin fragments (blacklist masked, >3000 frags) ---"
zcat Data/atac_filtered/${CARE}_fragments.tsv.gz | awk -F'\t' \
  -v CELLS=Data/cell_lists/${CARE}_cells.tsv \
  -v BLACKLIST=Data/ref/hg38-blacklist.v2.bed \
  -v W=500000 -v MINFRAGS=3000 -f Code/prep/bin_fragments.awk \
  2> Data/mintcna_input/logs/${CARE}_bin.log | gzip > Data/mintcna_input/${CARE}_counts_500kb.tsv.gz
cat Data/mintcna_input/logs/${CARE}_bin.log

echo "--- 2. ATAC + RNA input, config ---"
Rscript Code/prep/make_mintcna_input.R "$CARE"
Rscript Code/prep/make_rna_input.R     "$CARE"
Rscript Code/mintcna/make_config.R     "$CARE"
CFG=$WORK/Code/mintcna/config_${CARE}_multi.yaml

echo "--- 3. preprocess (ATAC + RNA l2r) ---"
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage preprocess

echo "--- 4. AG-DCAE ---"
OUT=$WORK/Result/MintCNA/${CARE}_multi
# Run the conda step in a SUBSHELL. `module load conda` conflicts with the R
# module and unloads it, and `conda deactivate` does not put it back -- so every
# Rscript after this point died with "command not found", after the 3 h of
# training had already succeeded. A subshell confines the module environment.
(
  module load conda
  source activate cae_cbam
  export PYTHONPATH="$REPO/python:${PYTHONPATH:-}"
  python "$REPO/scripts/run_denoise.py" \
    --input "$OUT/${CARE}_l2r_dna.csv" --output "$OUT/${CARE}_ae_dna.csv" \
    --latent "$OUT/${CARE}_latent_z.csv" --loss_plot "$OUT/${CARE}_loss.pdf" --config "$CFG"
)
command -v Rscript >/dev/null || { echo "R vanished from PATH after the conda step"; exit 1; }

echo "--- 5. align + scan ---"
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage align
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage scan

echo "--- 6. call, the configurations compared on P111T3 ---"
for spec in "0.10 0.20 a010e020" "0.10 0.25 a010e025" "0.05 0.30 a005e030"; do
  set -- $spec
  Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call \
    --alpha "$1" --st 3 --eps "$2" --tag "$3"
  cp "$OUT/${CARE}_mintcna_results.RData" "$OUT/${CARE}_results_$3.RData"
done
echo "=== $CARE done ==="
