#!/bin/bash
#SBATCH --job-name=crc_call_multi_d250
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=<REPRODUCE_DIR>/02_mintcna_runs/crc01/logs/call_multi_%j.out
#SBATCH --error=<REPRODUCE_DIR>/02_mintcna_runs/crc01/logs/call_multi_%j.err
set -uo pipefail; module load R/4.4; export R_LIBS_USER=<R_LIBS_USER>
cd <MINTCNA_ROOT>; t0=$(date +%s)
Rscript <REPRODUCE_DIR>/02_mintcna_runs/crc01/03_call_multi_d250.R > <REPRODUCE_DIR>/02_mintcna_runs/crc01/logs/call_multi_R.log 2>&1; rc=$?
echo "calling multi rc=$rc wall $(( $(date +%s) - t0 )) s  $(date)"; grep -E 'RUNTIME|saved' <REPRODUCE_DIR>/02_mintcna_runs/crc01/logs/call_multi_R.log; exit $rc
