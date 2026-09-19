#!/bin/bash
# Run from src/

tail -n +2 ../data/samplesheet.csv | cut -d, -f1 | while read -r s; do
    echo "=== $s ==="

    if [ ! -s ../data/raw/${s}_trimmed.fq ]; then
        echo "SKIP $s - no trimmed file" >&2
        continue
    fi

    rm -rf ../results/alignment/star/${s}_*

    STAR --genomeDir ../scerevisiae_star/ \
         --readFilesIn ../data/raw/${s}_trimmed.fq \
         --outFileNamePrefix ../results/alignment/star/${s}_ \
         --quantMode GeneCounts \
         --outSAMtype None \
         --runThreadN 4

    if [ ! -s ../results/alignment/star/${s}_ReadsPerGene.out.tab ]; then
        echo "FAILED: $s" >&2
    fi
done
