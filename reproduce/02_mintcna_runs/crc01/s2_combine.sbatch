#!/bin/bash
#SBATCH --job-name=crc_combine_d250
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=<REPRODUCE_DIR>/02_mintcna_runs/crc01/logs/combine_%j.out
#SBATCH --error=<REPRODUCE_DIR>/02_mintcna_runs/crc01/logs/combine_%j.err
set -uo pipefail; module load R/4.4; export R_LIBS_USER=<R_LIBS_USER>
cd <MINTCNA_ROOT> && Rscript <REPRODUCE_DIR>/02_mintcna_runs/crc01/02_combine_d250.R; echo "combine rc=$? $(date)"
