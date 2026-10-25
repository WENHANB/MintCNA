#!/bin/bash
# MintCNA on one CARE glioma 10x Multiome sample, ATAC + RNA and the ATAC-only ablation, exactly as run
# for the 20 samples of the paper: bin the fragments -> build the input -> preprocess -> AG-DCAE ->
# align -> change-point scan -> call (tag h15a005e020 = h 15, alpha 0.05, eps 0.20, 3 states).
#   bash reproduce/analysis/multiome/run_mintcna.sh <care_id>          (run from the repository root,
#                                                                       inside the mintcna conda env)
# Expects, under <MINTCNA_DATA>/ATAC_RNA (download_geo.sh ... make_rna_input.R in reproduce/analysis/multiome):
#   Data/atac_filtered/<care_id>_fragments.tsv.gz    fragments restricted to the authors' QC-passing nuclei
#   Data/cell_lists/<care_id>_cells.tsv              barcode, cell_type, nFrags, TSSEnrichment, is_reference
#   Data/ref/hg38-blacklist.v2.bed                   ENCODE blacklist
#   Data/gex_matrices/, Data/meta/sample_manifest.csv  for the RNA arm (make_rna_input.R)
set -euo pipefail
: "${MINTCNA_DATA:?set MINTCNA_DATA (see reproduce/README.md)}"
CARE=${1:?usage: run_multiome.sh <care_id>}
REPO=$(pwd)
CFG_DIR=$REPO/reproduce/analysis/multiome
WORK=$MINTCNA_DATA/ATAC_RNA
TAG=h15a005e020
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

cd "$WORK"
mkdir -p Data/mintcna_input/logs "Result/MintCNA/${CARE}_multi" "Result/MintCNA/${CARE}_atac"

echo "--- 1. bin fragments: 500 kb, blacklist masked, nuclei with > 3,000 fragments ---"
if [ ! -s "Data/mintcna_input/${CARE}_counts_500kb.tsv.gz" ]; then
  zcat "Data/atac_filtered/${CARE}_fragments.tsv.gz" | awk -F'\t' \
    -v CELLS="Data/cell_lists/${CARE}_cells.tsv" -v BLACKLIST=Data/ref/hg38-blacklist.v2.bed \
    -v W=500000 -v MINFRAGS=3000 -f "$REPO/scripts/atac_bin_fragments.awk" \
    2> "Data/mintcna_input/logs/${CARE}_bin.log" | gzip > "Data/mintcna_input/${CARE}_counts_500kb.tsv.gz"
fi

echo "--- 2. MintCNA input: ATAC arm (counts, bin annotation, cell table) and RNA arm ---"
Rscript "$REPO/scripts/make_atac_input.R" --counts "Data/mintcna_input/${CARE}_counts_500kb.tsv.gz" \
  --cells "Data/cell_lists/${CARE}_cells.tsv" --out "Data/mintcna_input/${CARE}" --genome hg38 --width 500000
Rscript "$REPO/reproduce/analysis/multiome/make_rna_input.R" "$CARE"

echo "--- 3. configs for this sample ---"
CFG=Result/MintCNA/config_${CARE}_multi.yaml;  sed "s/P111T2/${CARE}/g" "$CFG_DIR/config_multi.yaml" > "$CFG"
CFGA=Result/MintCNA/config_${CARE}_atac.yaml;  sed "s/P111T2/${CARE}/g" "$CFG_DIR/config_atac.yaml"  > "$CFGA"

echo "--- 4. ATAC + RNA: preprocess -> AG-DCAE -> align -> scan -> call ---"
OUT=Result/MintCNA/${CARE}_multi
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage preprocess --repo "$REPO"
PYTHONPATH="$REPO/MintCNA:${PYTHONPATH:-}" python "$REPO/scripts/run_denoise.py" --config "$CFG" \
  --input "$OUT/${CARE}_l2r_dna.csv" --output "$OUT/${CARE}_ae_dna.csv" \
  --latent "$OUT/${CARE}_latent_z.csv" --loss_plot "$OUT/${CARE}_loss.pdf"
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage align --repo "$REPO"
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage scan  --repo "$REPO"
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFG" --stage call  --repo "$REPO" --tag "$TAG"
cp "$OUT/${CARE}_mintcna_results.RData" "$OUT/${CARE}_results_${TAG}.RData"

echo "--- 5. ATAC only (ablation): same AG-DCAE output, RNA arm removed ---"
OUTA=Result/MintCNA/${CARE}_atac
for n in ae_dna.csv ref_dna.csv l2r_dna.csv latent_z.csv; do
  ln -sfn "$(readlink -f "$OUT/${CARE}_$n")" "$OUTA/${CARE}_$n"
done
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFGA" --stage align --repo "$REPO"
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFGA" --stage scan  --repo "$REPO"
Rscript "$REPO/scripts/run_mintcna.R" --config "$CFGA" --stage call  --repo "$REPO" --tag "$TAG"
cp "$OUTA/${CARE}_mintcna_results.RData" "$OUTA/${CARE}_results_${TAG}.RData"
echo "=== $CARE done: $OUT/${CARE}_results_${TAG}.RData, $OUTA/${CARE}_results_${TAG}.RData ==="
