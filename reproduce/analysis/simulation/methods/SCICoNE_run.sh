#!/bin/bash
# SCICoNE 1.0.0 on one simulation setting (github.com/cbg-ethz/SCICoNE).
#   bash reproduce/analysis/simulation/methods/SCICoNE_run.sh <setting>        (run from the repository root)
# Needs SCICONE_BIN (directory holding the breakpoint_detection and inference binaries and the python of
# the SCICoNE environment) and SCICONE_REPO (a clone of the SCICoNE repository, for scripts/segment_counts.py).
set -euo pipefail
: "${MINTCNA_DATA:?set MINTCNA_DATA (see reproduce/README.md)}"
: "${SCICONE_BIN:?set SCICONE_BIN}"; : "${SCICONE_REPO:?set SCICONE_REPO}"
BASE=${1:-c3_N100_100kb_l1e4}
SCIDIR=$MINTCNA_DATA/Simulation/cnasim/$BASE/Calling/SCICoNE
INDIR=$SCIDIR/input
mkdir -p "$INDIR" "$SCIDIR/logs"
PY=$SCICONE_BIN/python

# 1. counts -> cells x bins matrix
Rscript reproduce/analysis/simulation/methods/SCICoNE_input.R "$BASE"
IN=$INDIR/Y_cells_by_bins.csv
read N_CELLS N_BINS < <($PY - <<PY
import numpy as np, pathlib
f = "$IN"; first = pathlib.Path(f).read_text().splitlines()[0]
x = np.loadtxt(f, delimiter="," if "," in first else None); print(x.shape[0], x.shape[1])
PY
)

# 2. breakpoint detection -> segmented counts -> tree inference (all SCICoNE defaults)
cd "$SCIDIR"
"$SCICONE_BIN/breakpoint_detection" --d_matrix_file="$IN" --n_bins="$N_BINS" --n_cells="$N_CELLS"
REGION_SIZES=$SCIDIR/_segmented_region_sizes.txt
$PY "$SCICONE_REPO/scripts/segment_counts.py" "$IN" "$REGION_SIZES"
SEG_COUNTS=$INDIR/Y_cells_by_bins_segmented_counts.txt
read _ N_REGIONS < <($PY - <<PY
import numpy as np, pathlib
f = "$SEG_COUNTS"; first = pathlib.Path(f).read_text().splitlines()[0]
x = np.loadtxt(f, delimiter="," if "," in first else None); print(x.shape[0], x.shape[1])
PY
)
"$SCICONE_BIN/inference" --d_matrix_file="$SEG_COUNTS" --region_sizes_file="$REGION_SIZES" \
  --n_regions="$N_REGIONS" --n_cells="$N_CELLS" --postfix=run1
cd - > /dev/null

# 3. run1_inferred_cnvs.csv -> bins x cells state matrix (Calling/SCICoNE.RData)
Rscript reproduce/analysis/simulation/methods/SCICoNE_output.R "$BASE"
