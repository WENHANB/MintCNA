#!/bin/bash
#SBATCH --job-name=percell
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --cpus-per-task=4
#SBATCH --mem=96gb
#SBATCH --time=72:00:00
#SBATCH --output=<MINTCNA_ROOT>/ATAC_RNA/Logs/misc/percell_%j.out
# Memory is a first guess: the largest sample is ~4,900 malignant nuclei and the
# loader holds all six methods' cells x 27,000-bin matrices at once, ~1 GB each.
# Read MaxRSS from sacct afterwards and right-size, as with the benchmark jobs.
set -u
export TMPDIR=$HOME/rtmp; mkdir -p "$TMPDIR"
module load R/4.4
cd <MINTCNA_ROOT>/ATAC_RNA
Rscript Code/eval/pearson_percell.R
