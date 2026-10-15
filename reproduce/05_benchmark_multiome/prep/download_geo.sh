#!/bin/bash
# Download CARE IDH-mutant glioma multiome data from GEO.
#   GSE327580 = snATAC fragments   GSE326221 = snRNA count matrices
# Only samples flagged usable=yes in Data/meta/sample_manifest.csv are fetched
# (multiome + bulk DNA + >=500 ATAC nuclei -> 34 samples / 18 patients).
#
# Usage:  bash download_geo.sh {gex|atac|all}
# Safe to re-run: curl -C - resumes partial files, complete files are skipped.

set -u -o pipefail

ROOT="<MINTCNA_ROOT>/ATAC_RNA"
MANIFEST="$ROOT/Data/meta/sample_manifest.csv"
ATAC_DIR="$ROOT/Data/atac_fragments"
GEX_DIR="$ROOT/Data/gex_matrices"
LOG="$ROOT/Data/download.log"

# Individual supplementary files live under the per-SAMPLE tree, not the series
# one -- the series suppl dir only holds GSE*_RAW.tar, and filelist.txt merely
# describes that tar's contents. Series-level URLs 404.
#   https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9625nnn/GSM9625460/suppl/<file>
GEO_SAMPLES="https://ftp.ncbi.nlm.nih.gov/geo/samples"

# gsm_url <GSM123456> -> its suppl directory URL
gsm_url() {
    local gsm="$1"
    echo "$GEO_SAMPLES/${gsm:0:$(( ${#gsm} - 3 ))}nnn/$gsm/suppl"
}

mkdir -p "$ATAC_DIR" "$GEX_DIR"
MODE="${1:-all}"

say() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# fetch <url> <destdir> <filename> <expected_bytes>
fetch() {
    # NB: bash expands the whole `local` line before assigning, so dest must be
    # a separate statement or "$dir" is still unset when it is referenced.
    local url="$1" dir="$2" name="$3" want="$4"
    local dest="$dir/$name"
    if [[ -f "$dest" ]]; then
        local have
        have=$(stat -c%s "$dest" 2>/dev/null || echo 0)
        if [[ "$have" == "$want" ]]; then
            return 0                      # already complete
        fi
        say "resume $name ($have / $want bytes)"
    fi
    curl -sS -L --fail --retry 5 --retry-delay 5 --retry-all-errors \
         --connect-timeout 30 -C - -o "$dest" "$url/$name"
    local rc=$? have
    have=$(stat -c%s "$dest" 2>/dev/null || echo 0)
    if [[ $rc -ne 0 || "$have" != "$want" ]]; then
        say "FAIL $name (rc=$rc, got $have want $want)"
        return 1
    fi
    return 0
}

# Emit "field1<TAB>field2..." for usable rows. Columns are read by name so the
# script survives manifest column reordering.
usable_rows() {
    python3 - "$MANIFEST" <<'PY'
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    if r.get('usable') != 'yes':
        continue
    print('\t'.join([r['care_id'], r['atac_file'], r['atac_bytes'],
                     r['gex_files'], r['gsm_atac'], r['gsm_gex']]))
PY
}

# GEX file sizes are looked up from the GEO filelist rather than the manifest,
# which stores only their sum.
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
            if fetch "$gex_base" "$GEX_DIR" "$g" "${GEX_SIZE[$g]:-0}"; then
                n_ok=$((n_ok+1))
            else
                n_fail=$((n_fail+1))
            fi
        done
        say "gex  $care  ($gsm_gex)"
    fi
    if [[ "$MODE" == "atac" || "$MODE" == "all" ]]; then
        atac_base="$(gsm_url "$gsm_atac")"
        if fetch "$atac_base" "$ATAC_DIR" "$atac_file" "$atac_bytes"; then
            n_ok=$((n_ok+1))
            say "atac $care  $atac_file  ($(numfmt --to=iec "$atac_bytes"))"
        else
            n_fail=$((n_fail+1))
        fi
        # the tabix index is small and unlisted in the manifest
        fetch "$atac_base" "$ATAC_DIR" "${atac_file}.tbi.gz" \
              "$(awk -F'\t' -v f="${atac_file}.tbi.gz" '$2==f{print $4}' "$ROOT/Data/meta/filelist_GSE327580.txt")" || true
    fi
done < <(usable_rows)

say "=== done: $n_ok ok, $n_fail failed ==="
du -sh "$ATAC_DIR" "$GEX_DIR" 2>/dev/null | tee -a "$LOG"
