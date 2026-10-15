#!/bin/bash
#SBATCH --job-name=mintcna_rna
#SBATCH --output=<MINTCNA_ROOT>/Simulation/cnasim/c3_N100_100kb_l1e4/RNA_benchmark/logs/mintcna_rnaonly_%x_%j.out
#SBATCH --error=<MINTCNA_ROOT>/Simulation/cnasim/c3_N100_100kb_l1e4/RNA_benchmark/logs/mintcna_rnaonly_%x_%j.err
#SBATCH --account=feifeixiao
#SBATCH --qos=feifeixiao-b
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=08:00:00

module load R/4.4
BASE=${1:?need BASE arg}
echo "BASE=$BASE"
Rscript <REPRODUCE_DIR>/02_mintcna_runs/MintCNA_rnaonly.R "$BASE"
echo "mintcna_rnaonly $BASE exit=$?"
