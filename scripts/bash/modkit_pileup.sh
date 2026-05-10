#!/bin/bash
# =============================================
# modkit_pileup.sh
# Genero los archivos BEDMethyl para Control y Treatment
# =============================================

set -e
source scripts/bash/config.sh

echo "=== Ejecutando modkit pileup para Control ==="
modkit pileup --cpg \
  --mod-thresholds m:${MOD_THRESHOLD} --mod-thresholds h:${MOD_THRESHOLD} \
  --threads ${THREADS} \
  --ref ${REF_GENOME} \
  ${DATA_DIR}/Hlat.C.ONT.modmapped.bam \
  ${BAM_DIR}/Hlat.C.cpg.bed \
  --log-filepath ${BAM_DIR}/modkit_C.log

echo "=== Ejecutando modkit pileup para Treatment ==="
modkit pileup --cpg \
  --mod-thresholds m:${MOD_THRESHOLD} --mod-thresholds h:${MOD_THRESHOLD} \
  --threads ${THREADS} \
  --ref ${REF_GENOME} \
  ${DATA_DIR}/Hlat.T.ONT.modmapped.bam \
  ${BAM_DIR}/Hlat.T.cpg.bed \
  --log-filepath ${BAM_DIR}/modkit_T.log

echo "modkit pileup completado."
