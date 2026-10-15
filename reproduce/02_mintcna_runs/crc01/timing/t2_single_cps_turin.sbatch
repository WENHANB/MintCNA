#!/bin/bash
#SBATCH --job-name=crc_single_cps_timing
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --partition=hpg-turin
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=<MINTCNA_ROOT>/Code/CRC/rerun_d250/timing/logs/single_cps_turin_%j.out
#SBATCH --error=<MINTCNA_ROOT>/Code/CRC/rerun_d250/timing/logs/single_cps_turin_%j.err
set -uo pipefail; module load R/4.4; export R_LIBS_USER=<R_LIBS_USER>
cd <MINTCNA_ROOT>; echo "start $(date) host $(hostname)"; lscpu | grep -m1 "Model name"; t0=$(date +%s)
Rscript Code/CRC/rerun_d250/timing/04t_single_cps_timing.R > Code/CRC/rerun_d250/timing/logs/single_cps_turin_R.log 2>&1; rc=$?
echo "single cps TIMING rc=$rc wall $(( $(date +%s) - t0 )) s  $(date)"; grep -E 'TIMING' Code/CRC/rerun_d250/timing/logs/single_cps_turin_R.log; exit $rc
