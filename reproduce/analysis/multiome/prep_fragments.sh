#!/bin/bash
# Restrict each sample's Cell Ranger ARC fragments file to the authors' QC-passing nuclei
# (export_cell_lists.R), dropping ambient barcodes before any method reads the file.
#   bash reproduce/analysis/multiome/prep_fragments.sh [sample_list_file]      (default: the 20 samples)
# Input : ATAC_RNA/Data/atac_fragments/<file from the manifest>, ATAC_RNA/Data/cell_lists/<care_id>_barcodes.txt
# Output: ATAC_RNA/Data/atac_filtered/<care_id>_fragments.tsv.gz
# Run with MINTCNA_DATA set (see reproduce/README.md).
set -u -o pipefail
: "${MINTCNA_DATA:?set MINTCNA_DATA (see reproduce/README.md)}"
ROOT="${MINTCNA_DATA}/ATAC_RNA"
OUT="$ROOT/Data/atac_filtered"; mkdir -p "$OUT"
LIST="P59T2 P59T3 P104T1 P104T2 P106T1 P106T2 P100T1 P100T3 P103T1 P103T2 P105T1 P105T3 P108T1 P108T2 P109T1 P109T2 P110T1 P110T2 P111T2 P111T3"
[[ -n "${1:-}" ]] && LIST=$(cat "$1")

for care in $LIST; do
    atac_file=$(python3 - "$ROOT" "$care" <<'PY'
import csv, sys
root, care = sys.argv[1], sys.argv[2]
for r in csv.DictReader(open(f"{root}/Data/meta/sample_manifest.csv")):
    if r['care_id'] == care:
        print(r['atac_file']); break
PY
    )
    src="$ROOT/Data/atac_fragments/$atac_file"
    bc="$ROOT/Data/cell_lists/${care}_barcodes.txt"
    dest="$OUT/${care}_fragments.tsv.gz"
    if [[ -s "$dest" ]]; then echo "skip $care (exists)"; continue; fi
    [[ -f "$src" && -f "$bc" ]] || { echo "MISSING input for $care"; continue; }
    # keep the fragment rows whose barcode (column 4) is on the whitelist; drop comment lines
    zcat "$src" | awk -v BC="$bc" 'BEGIN{while((getline l < BC)>0) keep[l]=1} !/^#/ && ($4 in keep)' | gzip -c > "$dest"
    echo "$care  cells=$(wc -l < "$bc")  fragments_kept=$(zcat "$dest" | wc -l)"
done
