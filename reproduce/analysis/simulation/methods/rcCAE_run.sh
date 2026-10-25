#!/bin/bash
# rcCAE on one simulation setting (github.com/zhyu-lab/rccae, clone of 2026-03; Python 3.8 env + MATLAB).
#   bash reproduce/analysis/simulation/methods/rcCAE_run.sh <setting>        (run from the repository root)
# Needs RCCAE_REPO (the rccae clone) and an activated rcCAE Python environment; MATLAB on PATH.
set -euo pipefail
: "${MINTCNA_DATA:?set MINTCNA_DATA (see reproduce/README.md)}"
: "${RCCAE_REPO:?set RCCAE_REPO to the rccae clone}"
BASE=${1:-c3_N100_100kb_l1e4}
K=6                                           # number of SCHMM states
DIR=$MINTCNA_DATA/Simulation/cnasim/$BASE/Calling/rcCAE
INFILE=$DIR/rcCAE_input.txt
OUTDIR=$DIR/out
mkdir -p "$OUTDIR"

# 1. counts -> rcCAE input (bins x cells with gc and mappability columns)
Rscript reproduce/analysis/simulation/methods/rcCAE_input.R "$BASE"

# 2. rcCAE: train the autoencoder, then the SCHMM in MATLAB
( cd "$RCCAE_REPO"
  python ./cae/train.py --input "$INFILE" --output "$OUTDIR" \
    --epochs 100 --batch_size 64 --lr 0.0001 --latent_dim 3 --seed 0
  matlab -nodisplay -nosplash -r "addpath(genpath('${RCCAE_REPO}/hmm')); SCHMM('${OUTDIR}/lrc.txt','${OUTDIR}',${K}); exit"
)

# 3. copynumber.csv + ploidy.csv -> bins x cells state matrix (Calling/rcCAE.RData)
Rscript reproduce/analysis/simulation/methods/rcCAE_output.R "$BASE"
