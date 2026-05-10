#!/bin/bash
# =============================================
# run_pipeline_core.sh
# Ejecuta todos los pasos bash del pipeline
# =============================================

set -e
source scripts/bash/config.sh

echo "=== Iniciando Pipeline Core ==="

mkdir -p ${BAM_DIR} ${ANNOTATION_DIR}

./scripts/bash/samtools_stats.sh
./scripts/bash/modkit_pileup.sh
./scripts/bash/filter_coverage.sh
./scripts/bash/common_sites.sh
./scripts/bash/prepare_dss_input.sh
./scripts/bash/generate_bedgraphs.sh
./scripts/bash/genome_sizes.sh

echo "=== Pipeline Core completado correctamente ==="
echo "Todas las salidas están en ${BAM_DIR}/"
