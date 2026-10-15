#!/bin/bash
#SBATCH --job-name=crc_ae_d250_timing
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --partition=hpg-turin
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=<MINTCNA_ROOT>/Code/CRC/rerun_d250/timing/logs/ae_turin_%j.out
#SBATCH --error=<MINTCNA_ROOT>/Code/CRC/rerun_d250/timing/logs/ae_turin_%j.err
set -uo pipefail; module load pytorch/2.8.0
echo "start $(date) host $(hostname) partition hpg-turin"; lscpu | grep -m1 "Model name"; t0=$(date +%s)
python <MINTCNA_ROOT>/Code/CRC/rerun_d250/timing/01t_ae_timing.py; rc=$?
echo "AE TIMING rc=$rc  $(( $(date +%s) - t0 )) s  $(date)"; exit $rc
