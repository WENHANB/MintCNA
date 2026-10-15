#!/bin/bash
#SBATCH --job-name=crc_ae_d250
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=<MINTCNA_ROOT>/Code/CRC/rerun_d250/logs/ae_%j.out
#SBATCH --error=<MINTCNA_ROOT>/Code/CRC/rerun_d250/logs/ae_%j.err
set -uo pipefail; module load pytorch/2.8.0
echo "start $(date) host $(hostname) gpus=${CUDA_VISIBLE_DEVICES:-none}"; t0=$(date +%s)
python <MINTCNA_ROOT>/Code/CRC/rerun_d250/01_ae_d250.py; rc=$?
echo "AE rc=$rc  $(( $(date +%s) - t0 )) s  $(date)"; ls -la <MINTCNA_ROOT>/Data/CRC/DNAmeth/processed_l2r/*_d250_197.csv; exit $rc
