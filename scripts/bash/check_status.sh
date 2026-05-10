#!/bin/bash
# =============================================
# check_status.sh
# Verifica qué pasos ya están hechos y qué falta
# =============================================

source scripts/bash/config.sh

echo "=== Estado actual del pipeline ==="

echo "1. Archivos BAM originales:"
ls -lh ${DATA_DIR}/*.bam 2>/dev/null || echo "   No encontrados"

echo -e "\n2. Archivos BEDMethyl (modkit):"
ls -lh ${BAM_DIR}/*.cpg.bed 2>/dev/null || echo "   No generados aún"

echo -e "\n3. Archivos filtrados (cov10):"
ls -lh ${BAM_DIR}/*.cov10.bed 2>/dev/null | head -8 || echo "   No generados"

echo -e "\n4. Archivos DSS preparados:"
ls -lh ${BAM_DIR}/dss.*.txt 2>/dev/null || echo "   No generados"

echo -e "\n5. BedGraphs para IGV:"
ls -lh ${BAM_DIR}/*.bedgraph 2>/dev/null || echo "   No generados"

echo -e "\nListo. Usa esto para saber por dónde vas."
