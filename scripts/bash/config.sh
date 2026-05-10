#!/bin/bash
# =============================================
# config.sh - Configuración central del pipeline
# =============================================

# Rutas principales
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$( cd "${SCRIPT_DIR}/../../" && pwd )"
ANALYSIS_DIR="${BASE_DIR}/analisis"
BAM_DIR="${ANALYSIS_DIR}/bam_files"
ANNOTATION_DIR="${ANALYSIS_DIR}/annotation"
REFERENCE_DIR="${BASE_DIR}/reference_genome"
DATA_DIR="${BASE_DIR}/datos_secuenciacion_ori"

# Archivos importantes
REF_GENOME="${REFERENCE_DIR}/Heliconius_erato_lativitta_v1_-_scaffolds.fa"
GFF_FILE="${REFERENCE_DIR}/Hlat.v1.1.annotation.CAT.gff3.gz"

# Parámetros comunes
THREADS=14
COVERAGE_MIN=10
MOD_THRESHOLD=0.79

echo "→ Configuración cargada correctamente desde config.sh"
