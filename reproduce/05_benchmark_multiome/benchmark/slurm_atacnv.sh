#!/bin/bash
#SBATCH --job-name=atacnv
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --array=1-20%10
#SBATCH --cpus-per-task=2
#SBATCH --mem=16gb
#SBATCH --time=72:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Logs/AtaCNV/atacnv_%A_%a.out
set -u
export TMPDIR=$HOME/rtmp; mkdir -p "$TMPDIR"
module load R/4.4
ROOT=<MINTCNA_ROOT>/ATAC_RNA
SAMPLES=(P59T2 P59T3 P104T1 P104T2 P106T1 P106T2 P100T1 P100T3 P103T1 P103T2 \
         P105T1 P105T3 P108T1 P108T2 P109T1 P109T2 P110T1 P110T2 P111T2 P111T3)
CARE=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
echo "task $SLURM_ARRAY_TASK_ID -> $CARE"
cd "$ROOT" && Rscript Code/benchmark/run_atacnv.R "$CARE"
