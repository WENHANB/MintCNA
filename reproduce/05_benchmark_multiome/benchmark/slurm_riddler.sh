#!/bin/bash
#SBATCH --job-name=riddler
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --array=1-20%8
#SBATCH --cpus-per-task=2
#SBATCH --mem=24gb
#SBATCH --time=72:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Logs/RIDDLER/riddler_%A_%a.out
# window_matrix() loads the whole fragments file as a data.frame and then makes a
# filtered copy, so memory is sized for the largest sample (954 MB gzipped).
set -u
export TMPDIR=$HOME/rtmp; mkdir -p "$TMPDIR"
module load R/4.4
ROOT=<MINTCNA_ROOT>/ATAC_RNA
SAMPLES=(P59T2 P59T3 P104T1 P104T2 P106T1 P106T2 P100T1 P100T3 P103T1 P103T2 \
         P105T1 P105T3 P108T1 P108T2 P109T1 P109T2 P110T1 P110T2 P111T2 P111T3)
CARE=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
RES=${RES:-1MB}
echo "task $SLURM_ARRAY_TASK_ID -> $CARE at $RES"
cd "$ROOT" && Rscript Code/benchmark/run_riddler.R "$CARE" "$RES"
