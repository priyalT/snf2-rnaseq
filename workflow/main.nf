#!/usr/bin/env nextflow

include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trimgalore.nf'
include { STAR_INDEX } from './modules/star_index.nf'
include { STAR_ALIGN } from './modules/star_align.nf'
include { COUNT_MATRIX } from './modules/count_matrix.nf'
include { DESEQ2 } from './modules/deseq2.nf'
include { MULTIQC } from './modules/multiqc.nf'


 params {
    input: Path
    star_index: Path
    genomereads: Path
    genomeannotation: Path
    samplesheet: Path
    report_id: String
    datadir: Path
 }


workflow {
    main:
    read_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> file("${params.datadir}/${row.sample}.fastq")}
    FASTQC(read_ch)

    TRIM_GALORE(read_ch)
    if (!params.star_index) {
    genome_ch = channel.fromPath(params.genomereads)
    gtf_ch = channel.fromPath(params.genomeannotation)
    STAR_INDEX(genome_ch, gtf_ch)}
    star_index_ch = params.star_index
        ? channel.fromPath(params.star_index).first()
        : STAR_INDEX.out.index.first()
    
    STAR_ALIGN(TRIM_GALORE.out.trimmed_reads, star_index_ch)
 
    multiqc_files_ch = channel.empty().mix(
        FASTQC.out.zip,
        FASTQC.out.report,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        STAR_ALIGN.out.align_log,
    )
    multiqc_files_list = multiqc_files_ch.collect()
    MULTIQC(multiqc_files_list, params.report_id)


    samplesheet_ch = channel.fromPath(params.input).first()
    COUNT_MATRIX(STAR_ALIGN.out.gene_counts.collect(), samplesheet_ch)
    DESEQ2(COUNT_MATRIX.out.count_matrix, samplesheet_ch)


    publish:
    fastqc_zip = FASTQC.out.zip
    fastqc_html = FASTQC.out.report
    trimmed_reads = TRIM_GALORE.out.trimmed_reads
    trimming_reports = TRIM_GALORE.out.trimming_reports
    trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    gene_counts = STAR_ALIGN.out.gene_counts
    align_log = STAR_ALIGN.out.align_log
    count_matrix = COUNT_MATRIX.out.count_matrix
    deseq2_results = DESEQ2.out.deseq2_output
    multiqc_report = MULTIQC.out.report
    multiqc_data = MULTIQC.out.data


}

output {
    fastqc_zip {
        path 'fastqc'
    }
    fastqc_html {
        path 'fastqc'
    }
    trimmed_reads {
        path 'trimming'
    }
    trimming_reports {
        path 'trimming'
    }
    trimming_fastqc {
        path 'trimming'
    }
    gene_counts {
        path 'alignment'
    }
    align_log {
        path 'alignment'    
    }
    count_matrix {
        path 'counts'
    }
    deseq2_results {
        path 'deseq2'
    }
    multiqc_report {
        path 'multiqc'
    }
    multiqc_data {
        path 'multiqc'
    }

}