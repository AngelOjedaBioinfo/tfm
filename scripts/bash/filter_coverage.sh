#!/bin/bash
# =============================================
# filter_coverage.sh
# Filtrado por cobertura mínima
#
# Usage:
#   bash filter_coverage.sh <input.bed> <output_prefix> <min_coverage>
#
# Example:
#   bash filter_coverage.sh analisis/bam_files/Hlat.T.ds643.cpg.bed \
#                           analisis/bam_files/Hlat.T.ds643 10
# =============================================
set -euo pipefail
source scripts/bash/config.sh

INPUT="${1:-${BAM_DIR}/Hlat.T.cpg.bed}"
PREFIX="${2:-${BAM_DIR}/Hlat.T}"
MIN_COV="${3:-$COVERAGE_MIN}"

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input file not found: $INPUT" >&2
    exit 1
fi

OUT_5MC="${PREFIX}.5mC.cov${MIN_COV}.bed"
OUT_5HMC="${PREFIX}.5hmC.cov${MIN_COV}.bed"

echo "Filtrando sitios con cobertura >= ${MIN_COV}x..."
echo "  Input:    $INPUT"
echo "  5mC out:  $OUT_5MC"
echo "  5hmC out: $OUT_5HMC"

awk -v cov="$MIN_COV" '$4 == "m" && $10 >= cov' "$INPUT" > "$OUT_5MC" &
awk -v cov="$MIN_COV" '$4 == "h" && $10 >= cov' "$INPUT" > "$OUT_5HMC" &

wait

echo "Filtrado completado."
echo "  5mC sitios:  $(wc -l < "$OUT_5MC")"
echo "  5hmC sitios: $(wc -l < "$OUT_5HMC")"
