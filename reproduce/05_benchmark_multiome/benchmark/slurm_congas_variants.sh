#!/bin/bash
#SBATCH --job-name=congasVar
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --array=1-18
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=72:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Logs/CONGASplus/var_%A_%a.out
# Configuration search on three informative samples, refitting from the saved
# input checkpoint so the fragment binning is not repeated.
#   P105T1  strong oligodendroglioma, truth 1p -0.72 / 19q -0.71, currently K=1 flat CN3
#   P110T1  strongest truth in the cohort, 1p -0.92 / 19q -0.94
#   P59T2   astrocytoma negative control, truth flat at 1p/19q, real chr8q/chr9p events
# Grid: normal_cells {FALSE, TRUE} x lambda {0.5, 0.2, 0.8}
#   lambda weights the two likelihoods: toward 0 favours ATAC, toward 1 favours RNA.
set -u
export TMPDIR=$HOME/rtmp; mkdir -p "$TMPDIR"
module load R/4.4
ROOT=<MINTCNA_ROOT>/ATAC_RNA
export R_LIBS_USER=$ROOT/Tools/congas_env/rlib:$HOME/R/x86_64-pc-linux-gnu-library/4.4
export RETICULATE_PYTHON=$ROOT/Tools/congas_env/bin/python
SAMPLES=(P105T1 P110T1 P59T2)
NORMALS=(FALSE TRUE)
LAMBDAS=(0.5 0.2 0.8)
i=$((SLURM_ARRAY_TASK_ID - 1))
S=${SAMPLES[$((i % 3))]}
j=$((i / 3))
N=${NORMALS[$((j % 2))]}
L=${LAMBDAS[$((j / 2))]}
TAG="nc${N}_l${L}"
echo "task $SLURM_ARRAY_TASK_ID -> $S normal_cells=$N lambda=$L tag=$TAG"
cd "$ROOT" && Rscript Code/benchmark/refit_congas.R "$S" "$N" "$L" "$TAG"
