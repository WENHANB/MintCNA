#!/bin/bash
# Replace the path placeholders in every script under reproduce/ with your own locations.
#   bash reproduce/set_root.sh /path/to/project-root /path/to/R-library
# <MINTCNA_ROOT> = directory that holds Simulation/, Data/, Result/, ATAC_RNA/, figure4paper/
# <R_LIBS_USER>  = R library with the packages listed in environment.lock.yml
set -euo pipefail
ROOT=${1:?project root}; RLIB=${2:-$ROOT/Rlib}
grep -rlE '<MINTCNA_ROOT>|<R_LIBS_USER>|<BLUE_ROOT>' "$(dirname "$0")" | xargs -r sed -i -e "s|<MINTCNA_ROOT>|$ROOT|g" -e "s|<R_LIBS_USER>|$RLIB|g" -e "s|<BLUE_ROOT>|$ROOT/blue|g"
echo "placeholders replaced: ROOT=$ROOT  R_LIBS_USER=$RLIB"
