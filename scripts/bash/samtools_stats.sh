#!/bin/bash
# =============================================
# samtools_stats.sh
# Estadísticas básicas de calidad de los BAM
# =============================================

set -e
source scripts/bash/config.sh

echo "Generando estadísticas samtools..."

samtools stats ${DATA_DIR}/Hlat.C.ONT.modmapped.bam > ${BAM_DIR}/samtools_stats_C.txt 2> ${BAM_DIR}/samtools_stats_C.err
samtools stats ${DATA_DIR}/Hlat.T.ONT.modmapped.bam > ${BAM_DIR}/samtools_stats_T.txt 2> ${BAM_DIR}/samtools_stats_T.err

echo "Estadísticas samtools completadas."
