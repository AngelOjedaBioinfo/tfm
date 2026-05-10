#!/bin/bash
# =============================================
# extract_unnamed_genes.sh
# Extraigo secuencias de los 990 genes sin anotación funcional
# Para hacer BLAST contra SwissProt
# =============================================

set -e

GFF="reference_genome/Hlat.v1.1.annotation.CAT.gff3.gz"
UNNAMED_LIST="analisis/annotation/unnamed_gene_ids.txt"
BED_OUT="analisis/annotation/unnamed_genes.bed"
FASTA_OUT="analisis/annotation/unnamed_genes.fasta"

echo "Extrayendo coordenadas de genes sin anotación..."

zcat $GFF | grep -v "^#" | awk '$3 == "gene"' | \
while read line; do
    id=$(echo "$line" | grep -o 'gene_id=[^;]*' | cut -d= -f2)
    if grep -qx "$id" $UNNAMED_LIST; then
        chr=$(echo "$line" | cut -f1)
        start=$(echo "$line" | cut -f4)
        end=$(echo "$line" | cut -f5)
        echo -e "$chr\t$((start-1))\t$end\t$id"
    fi
done > $BED_OUT

echo "Coordenadas extraídas. Generando FASTA..."

bedtools getfasta \
  -fi reference_genome/Heliconius_erato_lativitta_v1_-_scaffolds.fa \
  -bed $BED_OUT \
  -name \
  -fo $FASTA_OUT

echo "Listo. Secuencias guardadas en $FASTA_OUT"
wc -l $BED_OUT
