#!/bin/bash
# Simulated scDNA-seq read counts: CNAsim 1.3.5 (conda env `cnasim`), 100 cells, 3 clones, 100 kb bins,
# at sequencing coverage 0.05x, 0.1x and 0.2x. CNAsim was run without a seed, so the generated read
# counts are not regenerated exactly; the nine deposited settings (Zenodo, see reproduce/README.md)
# are the ones used in the paper. The noise settings (_n0.05, _n0.1) are derived from these three
# by DNA_process.R <setting> <sigma>.
#   bash reproduce/analysis/simulation/cnasim.sh
set -euo pipefail
: "${MINTCNA_DATA:?set MINTCNA_DATA (see reproduce/README.md)}"
OUTROOT=$MINTCNA_DATA/Simulation/cnasim

for COV in 0.05 0.1 0.2; do
  case $COV in 0.1) OUT=$OUTROOT/c3_N100_100kb_l1e4 ;; *) OUT=$OUTROOT/c3_N100_100kb_l1e4_c$COV ;; esac
  mkdir -p "$OUT"
  cnasim -m 1 -o "$OUT" \
    -n 100 -n1 0.10 \
    -c 3 \
    -U \
    -p1 6 -k 10000 \
    -l1 10000 -l 5000000 \
    -B 100000 -C "$COV"
done
