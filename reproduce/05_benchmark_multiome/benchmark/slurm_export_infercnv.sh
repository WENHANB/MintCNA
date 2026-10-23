#!/bin/bash
#SBATCH --job-name=icnv_export
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=Log/icnv_export_%j.out
set -euo pipefail
cd <MINTCNA_ROOT>/ATAC_RNA
# Both configurations, so the headline (default) and the alternative (cut01) come
# from the same exporter. cut01 is re-exported because the script now reads
# run.final.infercnv_obj instead of the 22_denoise.* stage; they are the same
# object for that run, and re-exporting proves it rather than assuming it.
for TAG in default cut01_denoiseT; do
  echo "################ $TAG"
  INFERCNV_TAG=$TAG Tools/infercnv_env/bin/Rscript Code/benchmark/export_infercnv.R
done
