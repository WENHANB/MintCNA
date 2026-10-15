#!/bin/bash
# Filter each sample's fragments file down to the GEX called-cell barcodes.
#
# The GEO fragments are raw cellranger-arc output: one sample carries ~66,000
# barcodes in its first 12 MB alone, nearly all ambient. epiAneufinder reads the
# whole file into memory with fread() before applying any cell filter, so
# pre-filtering is what keeps 20 samples inside a sane memory footprint -- and it
# also removes ambient barcodes from every downstream tool at one stroke.
#
# Usage: bash prep_fragments.sh [sample_list_file]
#        default sample list = the 20 selected samples

set -u -o pipefail
ROOT="<MINTCNA_ROOT>/ATAC_RNA"
OUT="$ROOT/Data/atac_filtered"
LOG="$ROOT/Data/prep_fragments.log"
mkdir -p "$OUT"

SAMPLES=${1:-}
if [[ -z "$SAMPLES" ]]; then
    LIST="P59T2 P59T3 P104T1 P104T2 P106T1 P106T2 P100T1 P100T3 P103T1 P103T2 \
          P105T1 P105T3 P108T1 P108T2 P109T1 P109T2 P110T1 P110T2 P111T2 P111T3"
else
    LIST=$(cat "$SAMPLES")
fi

say() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
say "=== prep start: $(echo $LIST | wc -w) samples ==="

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
    # Whitelist is the AUTHORS' ArchR-QC'd cell list (doublets removed), not the
    # GEX called-cell list -- see export_cell_lists.R for why.
    bc="$ROOT/Data/cell_lists/${care}_barcodes.txt"
    dest="$OUT/${care}_fragments.tsv.gz"

    if [[ -s "$dest" ]]; then say "skip $care (exists)"; continue; fi
    [[ -f "$src" && -f "$bc" ]] || { say "MISSING input for $care"; continue; }

    n_bc=$(wc -l < "$bc")
    # Load the whitelist into an awk hash, then keep only fragment rows whose
    # barcode is in it. Comment lines are dropped: epiAneufinder names columns
    # positionally and does not skip '#'.
    zcat "$src" | awk -v BC="$bc" 'BEGIN{while((getline l < BC)>0) keep[l]=1}
        !/^#/ && ($4 in keep)' | gzip -c > "$dest"

    n_out=$(zcat "$dest" | wc -l)
    say "$care  cells=$n_bc  fragments_kept=$n_out  ($(awk -v a=$n_out -v b=$n_bc 'BEGIN{printf "%.0f", a/b}') per cell)  $(du -h "$dest" | cut -f1)"
done

say "=== prep done ==="
du -sh "$OUT"
