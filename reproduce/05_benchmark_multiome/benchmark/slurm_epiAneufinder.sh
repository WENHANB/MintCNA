#!/bin/bash
#SBATCH --job-name=epiAneu
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --array=41-80%12
#SBATCH --cpus-per-task=4
#SBATCH --mem=128gb
#SBATCH --time=72:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Logs/epiAneufinder/epiAneu_%A_%a.out

# 80 tasks = 20 samples x 4 configurations (windowSize x minFrags).
#
#   1-20  windowSize 1e5, minFrags 20000  -> w100kb_mf20000   DONE (job 39335320)
#  21-40  windowSize 1e5, minFrags  5000  -> w100kb_mf5000    DONE (job 39335320)
#  41-60  windowSize 1e6, minFrags 20000  -> w1mb_mf20000
#  61-80  windowSize 1e6, minFrags  5000  -> w1mb_mf5000
#
# windowSize is swept because the package gives conflicting guidance: it has no
# default and every doc passes 1e5, but the roxygen on the parameter itself says
# "Reccomended for sparse data - 1e6". At 100 kb this cohort has lambda = 0.45
# fragments per bin per nucleus, sparse by any definition, so 1e5 alone would test
# epiAneufinder in a regime its own authors flag as unsuitable.
#
# minFrags is swept because 20000 was tuned for standalone scATAC. Multiome splits
# the sequencing budget across two libraries from the same nuclei, so that gate
# keeps 23.4% of this cohort and drops 4 of 20 samples below 300 nuclei.
#
# Set --array=1-80%12 to run everything from scratch.

set -u
export TMPDIR=$HOME/rtmp
mkdir -p "$TMPDIR"
module load R/4.4

ROOT=<MINTCNA_ROOT>/ATAC_RNA
SAMPLES=(P59T2 P59T3 P104T1 P104T2 P106T1 P106T2 P100T1 P100T3 P103T1 P103T2 \
         P105T1 P105T3 P108T1 P108T2 P109T1 P109T2 P110T1 P110T2 P111T2 P111T3)

i=$((SLURM_ARRAY_TASK_ID - 1))
cfg=$((i / 20))
CARE=${SAMPLES[$((i % 20))]}

case $cfg in
  0) WINDOW=100000;  MINFRAGS=20000; TAG=w100kb_mf20000 ;;
  1) WINDOW=100000;  MINFRAGS=5000;  TAG=w100kb_mf5000  ;;
  2) WINDOW=1000000; MINFRAGS=20000; TAG=w1mb_mf20000   ;;
  3) WINDOW=1000000; MINFRAGS=5000;  TAG=w1mb_mf5000    ;;
  *) echo "bad config index $cfg"; exit 1 ;;
esac

echo "task $SLURM_ARRAY_TASK_ID -> $CARE  windowSize=$WINDOW  minFrags=$MINFRAGS  tag=$TAG"
cd "$ROOT"
Rscript Code/benchmark/run_epiAneufinder.R "$CARE" "$WINDOW" "$MINFRAGS" "$TAG"
