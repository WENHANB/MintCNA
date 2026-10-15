#!/bin/bash
#SBATCH --job-name=congasDiag
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --array=1-3
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=72:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Logs/CONGASplus/truthseg_%A_%a.out
# DIAGNOSTIC ONLY. Feeds CONGAS+ a truth-derived segmentation, which is circular
# and can never be reported as a benchmark score -- it exists solely to decide
# whether the arm-level segmentation is what is starving the model.
set -u
export TMPDIR=$HOME/rtmp; mkdir -p "$TMPDIR"
module load R/4.4
ROOT=<MINTCNA_ROOT>/ATAC_RNA
export R_LIBS_USER=$ROOT/Tools/congas_env/rlib:$HOME/R/x86_64-pc-linux-gnu-library/4.4
export RETICULATE_PYTHON=$ROOT/Tools/congas_env/bin/python
SAMPLES=(P105T1 P110T1 P59T2)
CARE=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
echo "DIAGNOSTIC: $CARE with truth-derived segmentation"
cd "$ROOT" && Rscript Code/benchmark/run_congas.R "$CARE" truth
