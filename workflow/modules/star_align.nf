#!/usr/bin/env nextflow

process STAR_ALIGN {

    container "community.wave.seqera.io/library/star:2.7.11b--5300af0cf0d14492"

    input:
    path trimmed_reads
    path star_index

    output:
    path "${trimmed_reads.simpleName}_ReadsPerGene.out.tab", emit: gene_counts
    path "${trimmed_reads.simpleName}_Log.final.out", emit: align_log


    script:
    """
    STAR --genomeDir ${star_index} \
         --readFilesIn ${trimmed_reads} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${trimmed_reads.simpleName}_ \
         --quantMode GeneCounts \
         --outSAMtype None \
         --runThreadN 4
    """
}