#!/bin/sh
#SBATCH --job-name=cnasim         # Job name
#SBATCH --output=log/cnasim%j.out
#SBATCH --error=log/cnasim%j.err
#SBATCH --account=feifeixiao	                  # your own sponser or account from this class
#SBATCH --qos=feifeixiao-b                     # your own sponser or account from this class
#SBATCH --mail-type=ALL                    # Mail events
#SBATCH --mail-user=<EMAIL>     # Where to send mail
#SBATCH --nodes=1                          # Use one node
#SBATCH --ntasks=1                      # Number of MPI ranks
#SBATCH --cpus-per-task=4               # Number of cores per MPI rank 
#SBATCH --mem=32gb                          # Total memory limit
#SBATCH --time=3-00:00:00                     # Time limit hrs:min:sec


module purge
module load conda
conda activate cnasim

BASE=<MINTCNA_ROOT>/Simulation
OUT=$BASE/cnasim/c3_N100_100kb_l1e4_c0.05
mkdir -p $OUT

cnasim -m 1 -o $OUT \
  -n 100 -n1 0.10 \
  -c 3 \
  -U \
  -p1 6 -k 10000\
  -l1 10000 -l 5000000\
  -B 100000 -C 0.05\