# snf2-rnaseq

A Nextflow-based RNA-seq analysis pipeline for differential gene expression, originally designed for comparing **snf2Δ** versus **wild type** in *Saccharomyces cerevisiae*.

## Overview

This pipeline processes bulk RNA-seq data end-to-end, from raw FASTQ reads to differential expression and functional enrichment.

![Nextflow Pipeline DAG](docs/dag.png)

The workflow consists of:

1. **Quality Control**: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. **Read Trimming**: [Trim Galore](https://github.com/FelixKrueger/TrimGalore)
3. **Alignment**: [STAR](https://github.com/alexdobin/STAR) (splice-aware alignment)
4. **Quantification**: STAR GeneCounts and custom count matrix aggregation
5. **Differential Expression**: [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
6. **Enrichment Analysis**: Functional enrichment of significant genes
7. **Reporting**: [MultiQC](https://multiqc.info/)

## Requirements

- [Nextflow](https://www.nextflow.io/) (>= 22.10.1)
- [Docker](https://www.docker.com/) (all dependencies are containerized)

## Usage

### Running the Test Profile

To verify the pipeline execution environment, run the included test profile:

```bash
nextflow run workflow/main.nf -profile test -with-dag docs/dag.png
```

### Full Pipeline Execution

To run the pipeline on your dataset, provide a sample sheet and reference genome files:

```bash
nextflow run workflow/main.nf \
    --input data/samplesheet.csv \
    --datadir data/raw \
    --genomereads genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
    --genomeannotation genome/Saccharomyces_cerevisiae.R64-1-1.114.gtf \
    --report_id snf2_analysis \
    -with-dag docs/dag.png \
    -resume
```

**Required Parameters:**
- `--input`: CSV sample sheet with a `sample` column.
- `--datadir`: Directory containing the raw `.fastq` files.
- `--genomereads`: Reference genome FASTA file.
- `--genomeannotation`: Reference genome GTF file.

## Repository Structure

- `workflow/`: Nextflow pipeline scripts, including `main.nf`, `modules/`, and `conf/`.
- `data/`: Sample sheets and raw data directory.
- `genome/`: Reference sequences and annotations.
- `report.qmd`: Quarto report template used for compiling downstream analysis and results.

## Output

Pipeline results are routed to specific output directories based on the process:
- `fastqc/` & `trimming/`: Pre- and post-trimming read quality reports.
- `alignment/`: STAR BAM files and alignment logs.
- `counts/`: Aggregated count matrix for all samples.
- `deseq2/`: Differential expression results (e.g., `de_snf2_vs_WT.csv`) and related plots.
- `enrichment/`: Functional enrichment results.
- `multiqc/`: A combined MultiQC HTML report summarizing the entire run.

## Background Data

The original analysis (documented in `report.qmd`) processed 12 samples (6 WT, 6 snf2Δ) from [ENA PRJEB5348](https://www.ebi.ac.uk/ena/browser/view/PRJEB5348). Snf2 is the catalytic ATPase of the SWI/SNF chromatin-remodelling complex, and its deletion perturbs transcription globally.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
