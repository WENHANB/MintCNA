#!/bin/bash
#SBATCH --job-name=congasSel
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --array=1-60%16
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=72:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Logs/CONGASplus/sel_%A_%a.out
# All 20 samples x normal_cells=TRUE x lambda {0.2,0.5,0.8}, WITH the package's
# own segments_selector_congas() step, which drops segments that are unimodal
# across cells. Omitting it was our error: it is the documented remedy for
# exactly the dilution we diagnosed (40 arms, most neutral).
# 8 cpus because segments_selector_congas fits one CONGAS+ model per segment and
# takes cores_ratio to parallelise them.
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
echo "task $SLURM_ARRAY_TASK_ID -> $S normal_cells=TRUE lambda=$L select=TRUE"
cd "$ROOT" && Rscript Code/benchmark/refit_congas.R "$S" TRUE "$L" "ncTRUE_sel_l${L}" TRUE
