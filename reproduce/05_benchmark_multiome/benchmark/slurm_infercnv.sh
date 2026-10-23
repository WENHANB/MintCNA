#!/bin/bash
#SBATCH --job-name=infercnv
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --array=1-20%8
#SBATCH --cpus-per-task=4
#SBATCH --mem=48gb
#SBATCH --time=72:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Logs/inferCNV/infercnv_%A_%a.out
# inferCNV lives in its own conda env: it needs rjags, which needs JAGS 4.x, and
# the cluster's jags module is an empty shell while conda-forge's default jags is
# 5.0.0 (rjags requires 4.x). bioconductor-infercnv from bioconda resolves the
# whole chain consistently.
set -u
export TMPDIR=$HOME/rtmp; mkdir -p "$TMPDIR"
ROOT=<MINTCNA_ROOT>/ATAC_RNA
RS=$ROOT/Tools/infercnv_env/bin/Rscript
SAMPLES=(P59T2 P59T3 P104T1 P104T2 P106T1 P106T2 P100T1 P100T3 P103T1 P103T2 \
         P105T1 P105T3 P108T1 P108T2 P109T1 P109T2 P110T1 P110T2 P111T2 P111T3)
CARE=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
echo "task $SLURM_ARRAY_TASK_ID -> $CARE"
cd "$ROOT" && "$RS" Code/benchmark/run_infercnv.R "$CARE"
