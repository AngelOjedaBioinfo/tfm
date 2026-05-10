#!/bin/bash
# =============================================
# common_sites.sh
# Intersección de sitios comunes entre Control y Treatment
# Solo los que tienen cobertura >=10x en ambas condiciones
# =============================================

set -e
source scripts/bash/config.sh

cd ${BAM_DIR}

echo "Calculando sitios comunes entre Control y Treatment..."

# 5mC
bedtools intersect -a Hlat.C.5mC.cov10.bed -b Hlat.T.5mC.cov10.bed -wa > common.5mC.C.bed
bedtools intersect -a Hlat.T.5mC.cov10.bed -b Hlat.C.5mC.cov10.bed -wa > common.5mC.T.bed

# 5hmC
bedtools intersect -a Hlat.C.5hmC.cov10.bed -b Hlat.T.5hmC.cov10.bed -wa > common.5hmC.C.bed
bedtools intersect -a Hlat.T.5hmC.cov10.bed -b Hlat.C.5hmC.cov10.bed -wa > common.5hmC.T.bed

echo "Sitios comunes calculados:"
wc -l common.*.bed
