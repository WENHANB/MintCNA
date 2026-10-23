#!/bin/bash
#SBATCH --job-name=congas
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --array=1-20%8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=72:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Logs/CONGASplus/congas_%A_%a.out
set -u
export TMPDIR=$HOME/rtmp; mkdir -p "$TMPDIR"
module load R/4.4
ROOT=<MINTCNA_ROOT>/ATAC_RNA
# Rcongas lives in its own library alongside the conda env that carries pyro/torch
export R_LIBS_USER=$ROOT/Tools/congas_env/rlib:$HOME/R/x86_64-pc-linux-gnu-library/4.4
export RETICULATE_PYTHON=$ROOT/Tools/congas_env/bin/python
SAMPLES=(P59T2 P59T3 P104T1 P104T2 P106T1 P106T2 P100T1 P100T3 P103T1 P103T2 \
         P105T1 P105T3 P108T1 P108T2 P109T1 P109T2 P110T1 P110T2 P111T2 P111T3)
CARE=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
echo "task $SLURM_ARRAY_TASK_ID -> $CARE"
cd "$ROOT" && Rscript Code/benchmark/run_congas.R "$CARE"
