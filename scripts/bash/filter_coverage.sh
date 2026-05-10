#!/bin/bash
# =============================================
# filter_coverage.sh
# Filtrado por cobertura mínima (>=10x)
# =============================================

set -e
source scripts/bash/config.sh

cd ${BAM_DIR}

echo "Filtrando sitios con cobertura >= ${COVERAGE_MIN}x..."

awk '$4 == "m" && $10 >= 10' Hlat.C.cpg.bed > Hlat.C.5mC.cov10.bed
awk '$4 == "h" && $10 >= 10' Hlat.C.cpg.bed > Hlat.C.5hmC.cov10.bed
awk '$4 == "m" && $10 >= 10' Hlat.T.cpg.bed > Hlat.T.5mC.cov10.bed
awk '$4 == "h" && $10 >= 10' Hlat.T.cpg.bed > Hlat.T.5hmC.cov10.bed

echo "Filtrado completado."
