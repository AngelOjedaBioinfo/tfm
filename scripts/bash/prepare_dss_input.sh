#!/bin/bash
# =============================================
# prepare_dss_input.sh
# Prepara los archivos de entrada para DSS (formato de 4 columnas)
# =============================================

set -e
source scripts/bash/config.sh

cd ${BAM_DIR}

echo "Preparando archivos de entrada para DSS..."

# 5mC
awk 'BEGIN{OFS="\t"} {print $1, $2+1, $10, $12}' common.5mC.C.bed > dss.5mC.C.txt
awk 'BEGIN{OFS="\t"} {print $1, $2+1, $10, $12}' common.5mC.T.bed > dss.5mC.T.txt

# 5hmC
awk 'BEGIN{OFS="\t"} {print $1, $2+1, $10, $12}' common.5hmC.C.bed > dss.5hmC.C.txt
awk 'BEGIN{OFS="\t"} {print $1, $2+1, $10, $12}' common.5hmC.T.bed > dss.5hmC.T.txt

echo "Archivos DSS preparados:"
ls -lh dss.*.txt
