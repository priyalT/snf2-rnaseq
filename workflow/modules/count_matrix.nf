process COUNT_MATRIX {

    container "community.wave.seqera.io/library/bioconductor-deseq2_r-tidyverse:958747a074ee937e"

    input:
    path gene_counts      
    path samplesheet      
    
    output:
    path "countmatrix.csv", emit: count_matrix

    script:
    """
    count_matrix.R ${samplesheet}
    """
}
