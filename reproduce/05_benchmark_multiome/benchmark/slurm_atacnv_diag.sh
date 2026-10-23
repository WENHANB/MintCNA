#!/bin/bash
#SBATCH --job-name=atacnvDiag
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --array=1-2
#SBATCH --cpus-per-task=2
#SBATCH --mem=16gb
#SBATCH --time=72:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Logs/AtaCNV/diag_%A_%a.out
set -u
export TMPDIR=$HOME/rtmp; mkdir -p "$TMPDIR"
module load R/4.4
ROOT=<MINTCNA_ROOT>/ATAC_RNA
SAMPLES=(P109T1 P111T3)
CARE=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
echo "diagnosing $CARE"
cd "$ROOT" && Rscript Code/benchmark/diagnose_atacnv_chr16.R "$CARE"
