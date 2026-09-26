#!/usr/bin/env nextflow

process FASTQC {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

    input: 
    path reads

    output: 
    path "${reads.simpleName}_fastqc.html", emit: report
    path "${reads.simpleName}_fastqc.zip", emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${reads}

    """
}