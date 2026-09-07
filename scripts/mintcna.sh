#!/usr/bin/env bash
# MintCNA - Top-level pipeline runner
#
# Orchestrates the full pipeline:
#   1. R: DNA (and optionally RNA) preprocessing
#   2. Python: AG-DCAE denoising of DNA L2R
#   3. R: align omics + segmentation + CNA calling
#
# Usage:
#   bash scripts/mintcna.sh --config config.yaml [--skip_denoise]

set -euo pipefail

# -------------------- Defaults --------------------
CONFIG=""
SKIP_DENOISE=0
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# -------------------- Parse args --------------------
while [[ $# -gt 0 ]]; do
  case $1 in
    -c|--config)
      CONFIG="$2"; shift 2 ;;
    --skip_denoise)
      SKIP_DENOISE=1; shift ;;
    -h|--help)
      grep '^#' "$0" | sed 's/^# \?//'
      exit 0 ;;
    *)
      echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$CONFIG" ]]; then
  echo "Error: --config <yaml> is required." >&2
  exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
  echo "Error: config file not found: $CONFIG" >&2
  exit 1
fi

# -------------------- Helpers --------------------
yaml_get() {
  # Minimal YAML extractor for simple key:value pairs (top-level or one nested level)
  local file="$1" key="$2"
  python3 -c "import yaml,sys; c=yaml.safe_load(open('$file')); ks='$key'.split('.'); v=c; \
[v := v.get(k) for k in ks]; print('' if v is None else v)" 2>/dev/null
}

OUT_DIR=$(yaml_get "$CONFIG" "output.dir")
PREFIX=$(yaml_get "$CONFIG" "output.prefix")
MODE=$(yaml_get "$CONFIG" "mode")

# Treat empty / None prefix as no prefix (no leading underscore)
if [[ -z "$PREFIX" || "$PREFIX" == "None" ]]; then
  P=""
else
  P="${PREFIX}_"
fi

# Auto-detect mode if blank
if [[ -z "$MODE" || "$MODE" == "None" ]]; then
  RNA=$(yaml_get "$CONFIG" "input.rna_counts")
  if [[ -n "$RNA" && "$RNA" != "None" ]]; then MODE="multi"; else MODE="dna"; fi
fi

mkdir -p "$OUT_DIR"
L2R_DNA="$OUT_DIR/${P}l2r_dna.csv"
AE_DNA="$OUT_DIR/${P}ae_dna.csv"
LATENT_Z="$OUT_DIR/${P}latent_z.csv"
LOSS_PLOT="$OUT_DIR/${P}loss_curve.pdf"

echo "============================================================"
echo "MintCNA Pipeline"
echo "  Config:    $CONFIG"
echo "  Repo:      $REPO_ROOT"
echo "  Mode:      $MODE"
echo "  Output:    $OUT_DIR"
echo "============================================================"

# -------------------- Stage 1: preprocessing --------------------
echo ""
echo ">>> [1/3] R preprocessing"
Rscript "$REPO_ROOT/scripts/run_mintcna.R" \
  --config "$CONFIG" --stage preprocess --repo "$REPO_ROOT"

# -------------------- Stage 2: AG-DCAE denoising --------------------
if [[ "$SKIP_DENOISE" -eq 1 ]]; then
  echo ""
  echo ">>> [2/3] Python denoising — SKIPPED (--skip_denoise)"
  cp "$L2R_DNA" "$AE_DNA"
else
  echo ""
  echo ">>> [2/3] Python AG-DCAE denoising"
  export PYTHONPATH="$REPO_ROOT/MintCNA:${PYTHONPATH:-}"
  python "$REPO_ROOT/scripts/run_denoise.py" \
    --input "$L2R_DNA" --output "$AE_DNA" \
    --latent "$LATENT_Z" --loss_plot "$LOSS_PLOT" \
    --config "$CONFIG"
fi

# -------------------- Stage 3: alignment + segmentation + calling --------------------
echo ""
echo ">>> [3/3] R alignment + segmentation + calling"
Rscript "$REPO_ROOT/scripts/run_mintcna.R" \
  --config "$CONFIG" --stage all_post_denoise --repo "$REPO_ROOT"

echo ""
echo "============================================================"
echo "MintCNA pipeline complete."
echo "Results: $OUT_DIR/${P}mintcna_results.RData"
echo "============================================================"
