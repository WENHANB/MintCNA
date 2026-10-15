#!/bin/bash
#SBATCH --job-name=congasNC
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --array=1-60%16
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=72:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Logs/CONGASplus/nc_%A_%a.out
# All 20 samples x normal_cells=TRUE x lambda {0.2, 0.5, 0.8}.
#
# The 3-sample search settled one thing and left one open. normal_cells=TRUE
# returns the median copy number to exactly 2.00 every time (2.57-2.94 without
# it), so it is adopted. Pearson improved from ~0 to ~0.1-0.2, and one sample
# reached 0.320 at lambda=0.8 -- but that was n=1, so lambda is swept across the
# full cohort rather than chosen from a single result.
set -u
export TMPDIR=$HOME/rtmp; mkdir -p "$TMPDIR"
module load R/4.4
ROOT=<MINTCNA_ROOT>/ATAC_RNA
export R_LIBS_USER=$ROOT/Tools/congas_env/rlib:$HOME/R/x86_64-pc-linux-gnu-library/4.4
export RETICULATE_PYTHON=$ROOT/Tools/congas_env/bin/python
SAMPLES=(P59T2 P59T3 P104T1 P104T2 P106T1 P106T2 P100T1 P100T3 P103T1 P103T2 \
         P105T1 P105T3 P108T1 P108T2 P109T1 P109T2 P110T1 P110T2 P111T2 P111T3)
LAMBDAS=(0.2 0.5 0.8)
i=$((SLURM_ARRAY_TASK_ID - 1))
S=${SAMPLES[$((i % 20))]}
L=${LAMBDAS[$((i / 20))]}
TAG="ncTRUE_l${L}"
echo "task $SLURM_ARRAY_TASK_ID -> $S normal_cells=TRUE lambda=$L"
cd "$ROOT" && Rscript Code/benchmark/refit_congas.R "$S" TRUE "$L" "$TAG"
