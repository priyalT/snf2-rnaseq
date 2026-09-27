#!/usr/bin/env nextflow

process STAR_INDEX {

    container "community.wave.seqera.io/library/star:2.7.11b--5300af0cf0d14492"

    input:
    path genomereads
    path genomeannotation

    output:
    path "star_index/", emit: index


    script:
    """
    mkdir -p star_index

    STAR --runMode genomeGenerate \
     --runThreadN ${task.cpus} \
     --genomeDir star_index \
     --genomeFastaFiles ${genomereads} \
     --sjdbGTFfile ${genomeannotation} \
     --sjdbOverhang 50 \
     --genomeSAindexNbases 11
    """
}