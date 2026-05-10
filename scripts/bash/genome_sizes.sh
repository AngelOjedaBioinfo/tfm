#!/bin/bash
# =============================================
# genome_sizes.sh
# Genera archivo con tamaños de cromosomas (para normalizar densidad de DMRs)
# =============================================

set -e
source scripts/bash/config.sh

cd ${BAM_DIR}

echo "Generando archivo de tamaños de cromosomas..."

awk '{print $1, $2}' ${REFERENCE_DIR}/Heliconius_erato_lativitta_v1_-_scaffolds.fa.fai | \
awk '$1 ~ /^Hel_chr/ {print $1, $2}' OFS='\t' > genome_sizes.txt

echo "Archivo generado: genome_sizes.txt"
head -n 8 genome_sizes.txt
