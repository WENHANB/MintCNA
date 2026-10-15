#!/bin/bash
#SBATCH --job-name=crc_call_single_d250
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=<MINTCNA_ROOT>/Code/CRC/rerun_d250/logs/call_single_%j.out
#SBATCH --error=<MINTCNA_ROOT>/Code/CRC/rerun_d250/logs/call_single_%j.err
set -uo pipefail; module load R/4.4; export R_LIBS_USER=<R_LIBS_USER>
cd <MINTCNA_ROOT>; t0=$(date +%s)
Rscript Code/CRC/rerun_d250/04_call_single_d250.R > Code/CRC/rerun_d250/logs/call_single_R.log 2>&1; rc=$?
echo "calling single rc=$rc wall $(( $(date +%s) - t0 )) s  $(date)"; grep -E 'RUNTIME|saved' Code/CRC/rerun_d250/logs/call_single_R.log; exit $rc
