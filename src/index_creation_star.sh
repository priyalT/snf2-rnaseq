#!/bin/bash
# Run from src/

mkdir -p ../scerevisiae_star

STAR --runMode genomeGenerate \
     --runThreadN 4 \
     --genomeDir ../scerevisiae_star \
     --genomeFastaFiles ../genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
     --sjdbGTFfile ../genome/Saccharomyces_cerevisiae.R64-1-1.114.gtf \
     --sjdbOverhang 50 \
     --genomeSAindexNbases 11
