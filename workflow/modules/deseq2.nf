process DESEQ2 {

    container "community.wave.seqera.io/library/bioconductor-apeglm_bioconductor-deseq2_pip_tidyverse:cee2c45c87144c2d"

    input:
    path count_matrix
    path samplesheet

    output:
    path "deseq2_output/", emit: deseq2_output

    script:
    """
    mkdir deseq2_output

    run_deseq2.R ${samplesheet}
    """

}