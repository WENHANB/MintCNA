#!/bin/bash
# MintCNA on CRC01 with the released package: DNA + RNA, then DNA only (the MintCNA^DNA arm of the paper).
#   bash reproduce/analysis/crc01/run_mintcna.sh
# Needs Data/CRC/mintcna_input/ (bin_counts.R, rna_counts.R), the mintcna conda environment, and runs
# from the root of the MintCNA repository with MINTCNA_DATA set.
# Output: Result/CRC/MintCNA/mintcna_results.RData and Result/CRC/MintCNA_dna/mintcna_results.RData
set -euo pipefail
: "${MINTCNA_DATA:?set MINTCNA_DATA (see reproduce/README.md)}"
REPO=$(pwd)
CFG=$REPO/reproduce/analysis/crc01/config_mintcna.yaml
cd "$MINTCNA_DATA"
[ -s Data/CRC/mintcna_input/DNA_RC.RData ] || { echo "missing Data/CRC/mintcna_input/DNA_RC.RData: run bin_counts.R first" >&2; exit 1; }
mkdir -p Result/CRC
sed -e 's#Result/CRC/MintCNA/#Result/CRC/MintCNA_dna/#' -e 's/^mode: "multi"/mode: "dna"/' \
    -e 's#^  rna_counts:.*#  rna_counts:      null#' -e 's#^  gene_annotation:.*#  gene_annotation: null#' \
    "$CFG" > Result/CRC/config_dna.yaml

bash "$REPO/scripts/mintcna.sh" --config "$CFG"
bash "$REPO/scripts/mintcna.sh" --config Result/CRC/config_dna.yaml
