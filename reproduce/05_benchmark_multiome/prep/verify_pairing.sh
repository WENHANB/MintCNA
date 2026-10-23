#!/bin/bash
# Verify same-nucleus pairing empirically for every downloaded sample.
#
# For each sample: take the GEX cell barcodes (the called-cell whitelist) and ask
# how many of them appear among the barcodes in the ATAC fragments file. 10x
# Multiome writes one shared barcode per nucleus, so a genuine multiome pair
# should overlap almost completely. Separate-aliquot assays would overlap only
# by chance on the shared 737K whitelist.
#
# Only the first N fragment lines are scanned (chr1 onward) -- enough for every
# real cell to appear at least once, and 100x cheaper than decompressing 52 GB.

set -u
ROOT="<MINTCNA_ROOT>/ATAC_RNA"
NLINES=${NLINES:-8000000}
OUT="$ROOT/Result/pairing_check.tsv"
mkdir -p "$ROOT/Result"

printf "care_id\tgex_barcodes\tatac_barcodes_scanned\toverlap\tpct\n" > "$OUT"

/apps/conda/25.7.0/bin/python3 - "$ROOT" <<'PY' | while IFS=$'\t' read -r care gexf atacf; do
import csv, sys
root=sys.argv[1]
for r in csv.DictReader(open(f"{root}/Data/meta/sample_manifest.csv")):
    if r['usable']!='yes': continue
    bc=[f for f in r['gex_files'].split(';') if f.endswith('barcodes.tsv.gz')][0]
    print('\t'.join([r['care_id'], bc, r['atac_file']]))
PY
    gex="$ROOT/Data/gex_matrices/$gexf"
    atac="$ROOT/Data/atac_fragments/$atacf"
    [[ -f "$gex" && -f "$atac" ]] || { echo "missing files for $care" >&2; continue; }

    read -r ngex natac nov < <(
        zcat "$gex" > /tmp/gex_$$.txt
        zcat "$atac" 2>/dev/null | head -n "$NLINES" | grep -v '^#' | cut -f4 | sort -u > /tmp/atac_$$.txt
        a=$(wc -l < /tmp/gex_$$.txt)
        b=$(wc -l < /tmp/atac_$$.txt)
        c=$(comm -12 <(sort /tmp/gex_$$.txt) /tmp/atac_$$.txt | wc -l)
        rm -f /tmp/gex_$$.txt /tmp/atac_$$.txt
        echo "$a $b $c"
    )
    pct=$(awk -v c="$nov" -v g="$ngex" 'BEGIN{printf "%.1f", (g?100*c/g:0)}')
    printf "%s\t%s\t%s\t%s\t%s\n" "$care" "$ngex" "$natac" "$nov" "$pct" >> "$OUT"
    echo "  $care  GEX=$ngex  overlap=$nov (${pct}%)"
done

echo "=> $OUT"
