#!/bin/bash
# Fill in the path placeholders in every script under reproduce/ and seed the data directories the figure
# scripts read from.
#   bash reproduce/set_root.sh /path/to/project-root [/path/to/R-library]
# <MINTCNA_ROOT>   = directory that holds Simulation/, Data/, Result/, ATAC_RNA/ and figure4paper/ (created here)
# <R_LIBS_USER>    = R library with the packages listed in environment.lock.yml (default: <root>/Rlib)
# <REPRODUCE_DIR>  = this folder (filled in automatically); the function libraries the scripts source are in lib/
set -euo pipefail
ROOT=$(cd "${1:?project root}" && pwd); RLIB=${2:-$ROOT/Rlib}; HERE=$(cd "$(dirname "$0")" && pwd)
grep -rlE '<MINTCNA_ROOT>|<R_LIBS_USER>|<BLUE_ROOT>|<REPRODUCE_DIR>' "$HERE" | xargs -r sed -i \
  -e "s|<MINTCNA_ROOT>|$ROOT|g" -e "s|<R_LIBS_USER>|$RLIB|g" -e "s|<BLUE_ROOT>|$ROOT/blue|g" -e "s|<REPRODUCE_DIR>|$HERE|g"
echo "placeholders replaced: ROOT=$ROOT  R_LIBS_USER=$RLIB  REPRODUCE_DIR=$HERE"
# The figure scripts read their source data from <root>/figure4paper/<figure>/data (or _work/data for the
# supplementary figures) and the merged simulation tables from <root>/Simulation/results/tables/merged.
# Seed those locations from data/ when they do not exist yet.
seed() { src=$HERE/data/$1; dst=$ROOT/$2; if [ -d "$src" ] && [ ! -e "$dst" ]; then mkdir -p "$dst"; cp -r "$src"/. "$dst"/; echo "seeded $dst"; fi; }
seed figure_source_data/figure2-radar          figure4paper/figure2-radar/data
seed figure_source_data/figure3-boundary-rmse  figure4paper/figure3-boundary-rmse/data
seed figure_source_data/figure4-CRC            figure4paper/figure4-CRC/data
seed figure_source_data/figure5-multiome       figure4paper/figure5-multiome/data
seed figure_source_data/supplementary          figure4paper/figureS-supplementary/_work/data
seed simulation_merged_tables                  Simulation/results/tables/merged
