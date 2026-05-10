#!/bin/bash
# =============================================
# generate_bedgraphs.sh
# Genera archivos bedGraph para visualización en IGV
# =============================================

set -e
source scripts/bash/config.sh

cd ${BAM_DIR}

echo "Generando bedGraphs para visualización..."

# 5mC
awk '$4 == "m" && $10 >= 10 {print $1"\t"$2"\t"$3"\t"$11}' Hlat.C.cpg.bed > Hlat.C.5mC.bedgraph
awk '$4 == "m" && $10 >= 10 {print $1"\t"$2"\t"$3"\t"$11}' Hlat.T.cpg.bed > Hlat.T.5mC.bedgraph

# 5hmC
awk '$4 == "h" && $10 >= 10 {print $1"\t"$2"\t"$3"\t"$11}' Hlat.C.cpg.bed > Hlat.C.5hmC.bedgraph
awk '$4 == "h" && $10 >= 10 {print $1"\t"$2"\t"$3"\t"$11}' Hlat.T.cpg.bed > Hlat.T.5hmC.bedgraph

echo "BedGraphs generados:"
ls -lh *.bedgraph
