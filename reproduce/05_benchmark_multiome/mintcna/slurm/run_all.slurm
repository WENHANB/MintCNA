#!/bin/bash
#SBATCH --job-name=mint_all
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --cpus-per-task=32
#SBATCH --mem=96gb
#SBATCH --output=<REPRODUCE_DIR>/05_benchmark_multiome/mintcna/slurm/mint_all_%A_%a.log

# The remaining 17 CARE samples at the LOCKED configuration. Nothing is swept
# here: h = 15, alpha = 0.10, eps = 0.20, st = 3, no genome smoothing, all read
# from config_<sample>_multi.yaml, which make_config.R copies from the P111T3
# template so a change cannot reach only some samples.
#
# Walltime is set per tier by the array's --time override at submission, because
# an over-long request cannot be backfilled: a 20 h ask sat in the queue for a
# day while the work itself takes one to five hours.
set -euo pipefail
SAMPLES=($SAMPLE_LIST)
CARE=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}

REPO=<REPRODUCE_DIR>/..
WORK=<MINTCNA_ROOT>/ATAC_RNA
cd "$WORK"
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=$WORK/Result/MintCNA/${CARE}_multi
CFG=$WORK/Code/mintcna/config_${CARE}_multi.yaml
echo "=== $CARE  ($(date)) ==="

if [ ! -s "Data/mintcna_input/${CARE}_counts_500kb.tsv.gz" ]; then
  echo "--- 1. bin fragments (blacklist masked, >3000 frags) ---"
  mkdir -p Data/mintcna_input/logs
  zcat Data/atac_filtered/${CARE}_fragments.tsv.gz | awk -F'\t' \
    -v CELLS=Data/cell_lists/${CARE}_cells.tsv \
    -v BLACKLIST=Data/ref/hg38-blacklist.v2.bed \
    -v W=500000 -v MINFRAGS=3000 -f Code/prep/bin_fragments.awk \
    2> Data/mintcna_input/logs/${CARE}_bin.log | gzip > Data/mintcna_input/${CARE}_counts_500kb.tsv.gz
  cat Data/mintcna_input/logs/${CARE}_bin.log
else echo "--- 1. counts already binned, skipping ---"; fi

echo "--- 2. ATAC + RNA input, config ---"
Rscript Code/prep/make_mintcna_input.R "$CARE"
Rscript Code/prep/make_rna_input.R     "$CARE"
Rscript <REPRODUCE_DIR>/05_benchmark_multiome/mintcna/make_config.R     "$CARE"

echo "--- 3. preprocess ---"
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage preprocess

echo "--- 4. AG-DCAE ---"
# Subshell: `module load conda` conflicts with the R module and unloads it, and
# `conda deactivate` does not restore it. Job 39636410 lost 3 h of finished
# training to "Rscript: command not found" on the next line.
(
  module load conda
  source activate cae_cbam
  export PYTHONPATH="$REPO/python:${PYTHONPATH:-}"
  python "$REPO/scripts/run_denoise.py" \
    --input "$OUT/${CARE}_l2r_dna.csv" --output "$OUT/${CARE}_ae_dna.csv" \
    --latent "$OUT/${CARE}_latent_z.csv" --loss_plot "$OUT/${CARE}_loss.pdf" --config "$CFG"
)
command -v Rscript >/dev/null || { echo "R vanished from PATH after the conda step"; exit 1; }

echo "--- 5. align + scan (checkpoint kept: ${CARE}_scan_h15.rds) ---"
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage align
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage scan

echo "--- 6. call ---"
TAG=h15a010e020
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call --tag "$TAG"
cp "$OUT/${CARE}_mintcna_results.RData" "$OUT/${CARE}_results_${TAG}.RData"
# the parallel filesystem produced one 1.1 MB file of nul bytes in an earlier run
Rscript -e "e<-new.env(); load('$OUT/${CARE}_results_${TAG}.RData', e);
            stopifnot(nrow(get('cna_state_matrix',e))>0); cat('verified\n')"
echo "=== $CARE done ($(date)) ==="
