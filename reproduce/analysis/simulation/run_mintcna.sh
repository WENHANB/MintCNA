#!/bin/bash
# MintCNA on one simulation setting: DNA + RNA, then DNA only (the MintCNA^DNA arm of the paper).
#   bash reproduce/analysis/simulation/run_mintcna.sh <setting>
# Needs <setting>/DNA_RC.RData and <setting>/mintcna_input/ (make_mintcna_input.R), the mintcna conda
# environment, and runs from the root of the MintCNA repository with MINTCNA_DATA set.
# Output: <setting>/mintcna_multi/mintcna_results.RData and <setting>/mintcna_dna/mintcna_results.RData
set -euo pipefail
: "${MINTCNA_DATA:?set MINTCNA_DATA (see reproduce/README.md)}"
SETTING=${1:?usage: run_mintcna.sh <setting>}
REPO=$(pwd)
TPL=$REPO/reproduce/analysis/simulation/config_mintcna.yaml
cd "$MINTCNA_DATA"
D=Simulation/cnasim/$SETTING
[ -s "$D/DNA_RC.RData" ] || { echo "missing $D/DNA_RC.RData: run make_mintcna_input.R first" >&2; exit 1; }

sed "s/__SETTING__/$SETTING/g" "$TPL" > "$D/config_multi.yaml"
sed -e 's#/mintcna_multi/#/mintcna_dna/#' -e 's/^mode: "multi"/mode: "dna"/' \
    -e 's#^  rna_counts:.*#  rna_counts:      null#' -e 's#^  gene_annotation:.*#  gene_annotation: null#' \
    "$D/config_multi.yaml" > "$D/config_dna.yaml"

bash "$REPO/scripts/mintcna.sh" --config "$D/config_multi.yaml"
bash "$REPO/scripts/mintcna.sh" --config "$D/config_dna.yaml"
