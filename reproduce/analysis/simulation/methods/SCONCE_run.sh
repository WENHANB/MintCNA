#!/bin/bash
# SCONCE on one simulation setting (github.com/NielsenBerkeleyLab/sconce, commit 1f1e1f9, built from source;
# needs gsl, boost and gcc).
#   bash reproduce/analysis/simulation/methods/SCONCE_run.sh <setting>        (run from the repository root)
# Needs SCONCE_REPO (the sconce clone, with bin/sconce built).
set -euo pipefail
: "${MINTCNA_DATA:?set MINTCNA_DATA (see reproduce/README.md)}"
: "${SCONCE_REPO:?set SCONCE_REPO to the sconce clone}"
BASE=${1:-c3_N100_100kb_l1e4}
MAXKPLOID=5
BASE_DIR=$MINTCNA_DATA/Simulation/cnasim/$BASE/Calling/SCONCE
SCONCE_BIN=$SCONCE_REPO/bin/sconce
mkdir -p "$BASE_DIR/ref" "$BASE_DIR/out" "$BASE_DIR/logs"

# 1. counts -> one BED of read depths per cell, plus the diploid / tumour file lists
Rscript reproduce/analysis/simulation/methods/SCONCE_input.R "$BASE"
DIP_LIST=$BASE_DIR/ref/diploidFileList.txt
TUM_LIST=$BASE_DIR/tumorFileList.txt
HEALTHY_AVG=$BASE_DIR/ref/healthy_avg.bed
MEANVAR=$BASE_DIR/ref/meanVar.txt

# 2. diploid reference and mean-variance relationship, then SCONCE on every tumour cell
Rscript "$SCONCE_REPO/scripts/avgDiploid.R"        "$DIP_LIST" "$HEALTHY_AVG"
Rscript "$SCONCE_REPO/scripts/fitMeanVarRlnshp.R"  "$DIP_LIST" "$MEANVAR"
while read -r tumorbed; do
  [ -z "$tumorbed" ] && continue
  cell=$(basename "$tumorbed" .bed)
  "$SCONCE_BIN" --diploid "$HEALTHY_AVG" --tumor "$tumorbed" --meanVarCoefFile "$MEANVAR" \
    --outputBase "$BASE_DIR/out/${cell}_k${MAXKPLOID}" --maxKploid "$MAXKPLOID" \
    > "$BASE_DIR/logs/${cell}_k${MAXKPLOID}.log" 2> "$BASE_DIR/logs/${cell}_k${MAXKPLOID}.err"
done < "$TUM_LIST"

# 3. per-cell BEDs -> bins x cells state matrix (Calling/SCONCE.RData)
Rscript reproduce/analysis/simulation/methods/SCONCE_output.R "$BASE"
