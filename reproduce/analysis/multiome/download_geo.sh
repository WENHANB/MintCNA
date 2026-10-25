#!/bin/bash
# Download the CARE IDH-mutant glioma 10x Multiome data from GEO: snATAC fragment files (GSE327580) and
# snRNA count matrices (GSE326221) of the samples flagged usable=yes in Data/meta/sample_manifest.csv.
#   bash reproduce/analysis/multiome/download_geo.sh {gex|atac|all}
# Output: ATAC_RNA/Data/atac_fragments/, ATAC_RNA/Data/gex_matrices/   (resumable; complete files are skipped)
# Expected file sizes come from the manifest and from the GEO file lists
# Data/meta/filelist_GSE326221.txt and filelist_GSE327580.txt. Run with MINTCNA_DATA set.
set -u -o pipefail
: "${MINTCNA_DATA:?set MINTCNA_DATA (see reproduce/README.md)}"
ROOT="${MINTCNA_DATA}/ATAC_RNA"
MANIFEST="$ROOT/Data/meta/sample_manifest.csv"
ATAC_DIR="$ROOT/Data/atac_fragments"
GEX_DIR="$ROOT/Data/gex_matrices"
LOG="$ROOT/Data/download.log"
GEO_SAMPLES="https://ftp.ncbi.nlm.nih.gov/geo/samples"      # per-sample supplementary files
MODE="${1:-all}"
mkdir -p "$ATAC_DIR" "$GEX_DIR"

say() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
gsm_url() { local gsm="$1"; echo "$GEO_SAMPLES/${gsm:0:$(( ${#gsm} - 3 ))}nnn/$gsm/suppl"; }

# fetch <url> <destdir> <filename> <expected_bytes>: resume a partial file, skip a complete one
fetch() {
    local url="$1" dir="$2" name="$3" want="$4"
    local dest="$dir/$name"
    if [[ -f "$dest" ]]; then
        local have; have=$(stat -c%s "$dest" 2>/dev/null || echo 0)
        [[ "$have" == "$want" ]] && return 0
        say "resume $name ($have / $want bytes)"
    fi
    curl -sS -L --fail --retry 5 --retry-delay 5 --retry-all-errors --connect-timeout 30 -C - -o "$dest" "$url/$name"
    local rc=$? have; have=$(stat -c%s "$dest" 2>/dev/null || echo 0)
    if [[ $rc -ne 0 || "$have" != "$want" ]]; then say "FAIL $name (rc=$rc, got $have want $want)"; return 1; fi
    return 0
}

# usable manifest rows: care_id, atac_file, atac_bytes, gex_files, gsm_atac, gsm_gex
usable_rows() {
    python3 - "$MANIFEST" <<'PY'
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    if r.get('usable') == 'yes':
        print('\t'.join([r['care_id'], r['atac_file'], r['atac_bytes'], r['gex_files'], r['gsm_atac'], r['gsm_gex']]))
PY
}

declare -A GEX_SIZE
while IFS=$'\t' read -r kind name _ size _; do
    [[ "$kind" == "File" ]] && GEX_SIZE["$name"]="$size"
done < "$ROOT/Data/meta/filelist_GSE326221.txt"

say "=== download start (mode=$MODE) ==="
n_ok=0; n_fail=0
while IFS=$'\t' read -r care atac_file atac_bytes gex_files gsm_atac gsm_gex; do
    if [[ "$MODE" == "gex" || "$MODE" == "all" ]]; then
        gex_base="$(gsm_url "$gsm_gex")"
        IFS=';' read -ra gfs <<< "$gex_files"
        for g in "${gfs[@]}"; do
            [[ -z "$g" ]] && continue
            if fetch "$gex_base" "$GEX_DIR" "$g" "${GEX_SIZE[$g]:-0}"; then n_ok=$((n_ok+1)); else n_fail=$((n_fail+1)); fi
        done
        say "gex  $care  ($gsm_gex)"
    fi
    if [[ "$MODE" == "atac" || "$MODE" == "all" ]]; then
        atac_base="$(gsm_url "$gsm_atac")"
        if fetch "$atac_base" "$ATAC_DIR" "$atac_file" "$atac_bytes"; then
            n_ok=$((n_ok+1)); say "atac $care  $atac_file  ($(numfmt --to=iec "$atac_bytes"))"
        else n_fail=$((n_fail+1)); fi
        fetch "$atac_base" "$ATAC_DIR" "${atac_file}.tbi.gz" \
              "$(awk -F'\t' -v f="${atac_file}.tbi.gz" '$2==f{print $4}' "$ROOT/Data/meta/filelist_GSE327580.txt")" || true
    fi
done < <(usable_rows)
say "=== done: $n_ok ok, $n_fail failed ==="
