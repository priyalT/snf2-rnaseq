process ENRICHMENT {

    container "community.wave.seqera.io/library/bioconductor-clusterprofiler_bioconductor-enrichplot_bioconductor-org.sc.sgd.db:014325c664b398e2"

    input:
    path deseq_dataset

    output:
    path "enrichment_output/", emit: enrichment_output

    script:
    """
    mkdir enrichment_output

    enrichment.R ${deseq_dataset}/de_snf2_vs_WT.csv
    """
}