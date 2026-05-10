#!/bin/bash
# =============================================
# run_blast.sh
# BLASTx de los genes sin anotación funcional (990 genes)
# Uso SwissProt porque tiene anotación curada
# =============================================

set -e

echo "Ejecutando blastx de genes sin anotación..."

blastx \
  -query analisis/annotation/unnamed_genes.fasta \
  -db blast_db/uniprot_sprot \
  -out analisis/annotation/unnamed_genes_blast.txt \
  -evalue 1e-5 \
  -num_threads 8 \
  -max_target_seqs 1 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle"

echo "BLAST terminado. Seleccionando mejor hit por gen..."

cd analisis/annotation
sort -k1,1 -k5,5g unnamed_genes_blast.txt | awk '!seen[$1]++' > blast_best_hits.txt

echo "Mejores hits guardados en blast_best_hits.txt"
