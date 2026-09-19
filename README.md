# snf2-rnaseq
## Summary

Bulk RNA-seq from raw reads: differential gene expression in *Saccharomyces cerevisiae* involving two conditions:
**snf2Δ** versus **wild type**.

Twelve samples were included (6 WT, 6 snf2Δ) from [ENA PRJEB5348](https://www.ebi.ac.uk/ena/browser/view/PRJEB5348)
taken end to end from: 

**raw FASTQ** → QC → adapter/quality trimming → splice-aware
alignment → gene-level counting → **differential expression**


## Scope

**In scope:** 12 of the 96 available samples, one contrast (snf2Δ vs WT), the
standard upstream workflow, with each parameter choice and its reasoning
recorded, and the known weaknesses listed rather than hidden.

**Not in scope:** novel biology, method development, production packaging, or the
replicate-depth question the source dataset was designed around.

## Result

The most significant change in the whole experiment is
**SNF2 itself**, the gene that was knocked out.

| Gene | Systematic ID | log2FC (shrunken) | baseMean | padj |
|---|---|---:|---:|---:|
| **SNF2** | YOR290C | **−5.35** | 120 | 1.5 × 10⁻⁵⁸ |

Summary of the full contrast (snf2Δ vs WT, Wald test, Benjamini–Hochberg):

| | |
|---|---:|
| Gene features in the annotation (Ensembl 114 on R64-1-1; 6,600 protein-coding) | 7,127 |
| Genes passing the count filter | 5,401 |
| Genes with a usable padj (30 dropped as Cook's-distance outliers) | 5,371 |
| Significant at padj < 0.05 | 1,348 |
| Significant **and** \|log2FC\| > 1 | 296 (73 up, 223 down) |

A quarter of the tested transcriptome moves at padj < 0.05. That is expected because
Snf2 is the catalytic ATPase of the SWI/SNF chromatin-remodelling complex, so
deleting it perturbs transcription globally rather than in one pathway. 

**Observed pattern:** Among the strongest downregulated genes
are `PHO84`, `PHO12`, `SPL2`, `VTC3` and `VTC4` which are all phosphate-responsive.
SWI/SNF is known to be required for chromatin remodelling at PHO promoters, so
this is a coherent signal rather than noise, and it was not something the
analysis was set up to look for.

### Figures

| | |
|---|---|
| `plots/PCA_plot.png` | VST-transformed samples; WT and snf2Δ separate on PC1 |
| `plots/MA_plot.png` | Unshrunken log2 fold changes |
| `plots/MA_shrunkenLFC_plot.png` | The same contrast after `apeglm` shrinkage |
| `plots/ShrunkenLFC_comparison_plot.png` | Unshrunken vs shrunken LFC, with y = x |
| `plots/Volcano_plot.png` | padj < 0.05 & \|log2FC\| > 1 highlighted |



## Data

[Gierliński et al. 2015](https://doi.org/10.1093/bioinformatics/btv425) /
[Schurch et al. 2016](https://doi.org/10.1261/rna.053959.115) generated 48
biological replicates each of wild-type and snf2Δ *S. cerevisiae*.

| | |
|---|---|
| Accession | ENA PRJEB5348 |
| Organism | *S. cerevisiae*, BY4741 background |
| Samples used | 12 of 96 (6 WT, 6 snf2Δ) — see `data/samplesheet.csv` |
| Reads | 51 bp, single-end, Illumina HiSeq 2000 (2014) |
| Reference | Ensembl R64-1-1, `dna.toplevel.fa` + release-114 GTF |
| Strandedness | **Unstranded** — verified: forward/reverse counts split 51.5% / 48.5% (475,612 vs 448,211), so column 2 of `ReadsPerGene.out.tab` is correct |

Only 12 of the 96 available samples were used because of storage constrictions (roughly 20 GB of free disk). Raw and trimmed FASTQ, the STAR
index and BAM-adjacent outputs are **not** committed (see `.gitignore`).

---

## Pipeline

```
data/raw/*.fastq
      │
      ├─ FastQC ──────────────► results/*_fastqc.html
      │                          │
      │                          └─ MultiQC ──► results/multiqc_report.html
      ▼
   Trim Galore  (Cutadapt + quality trim)
      │
      ▼  data/raw/*_trimmed.fq
   STAR --quantMode GeneCounts        (index: scerevisiae_star/)
      │
      ▼  results/alignment/star/*_ReadsPerGene.out.tab
   count_matrix.R    ── merge per-sample counts, drop N_* rows
      │
      ▼  results/countmatrix.csv
   run_deseq2.R      ── DESeq2: size factors → dispersions → Wald → BH → apeglm
      │
      ▼  results/de_snf2_vs_WT.csv  +  plots/
```

---

## Repository layout

```
├── data/
│   ├── ena_report.csv           # ENA metadata for the run accessions
│   ├── samplesheet.csv          # sample,condition
│   └── raw/                     # FASTQ (gitignored)
├── genome/                      # reference FASTA + GTF (gitignored)
├── scerevisiae_star/            # STAR index (gitignored)
├── src/
│   ├── fastqc.sh                # FastQC over raw reads
│   ├── trim.sh                  # Trim Galore
│   ├── index_creation_star.sh   # STAR genomeGenerate
│   ├── alignment.sh             # per-sample STAR, idempotent
│   ├── bam_exploration.sh       # ad-hoc SAM/BAM inspection
│   ├── count_matrix.R           # per-sample counts → count matrix
│   ├── run_deseq2.R             # normalisation, DE, shrinkage, plots
│   ├── explore.py               # FASTQ/quality-score exploration
│   └── explore.ipynb
├── results/
│   ├── countmatrix.csv          # 7,127 genes × 12 samples, raw integer counts
│   ├── de_snf2_vs_WT.csv        # full DE table, ordered by padj
│   └── multiqc_data*/
└── plots/
```

---

## Reproducing

### Requirements

These are the exact versions that produced the results in this repository, and
how each one was installed on the machine that ran it (macOS, Apple Silicon).


| Tool | Version | Installed via |
|---|---|---|
| STAR | 2.7.11b | compiled from the GitHub release |
| Trim Galore | 2.3.0 | Homebrew (`brew install trim-galore`) |
| FastQC | 0.12.1 | official zip, unpacked locally (gitignored) |
| MultiQC | 1.35 | uv, from `pyproject.toml` + `uv.lock` — `uv sync` |
| R | 4.6.0 | CRAN installer |
| DESeq2 | 1.52.0 | Bioconductor — `Rscript src/install_r_packages.R` |
| apeglm | 1.34.0 | Bioconductor — same script |

Trim Galore 2.x is a self-contained binary and no longer needs a separate
`cutadapt` install; 0.6.x did.

Pinning these into one reproducible environment is the job of a container, and
is deferred to the Nextflow work in [Future improvements](#future-improvements).
```bash
# 1. reference — Ensembl release 114, S. cerevisiae R64-1-1
#    unsoftmasked toplevel FASTA + matching GTF, both gunzipped
mkdir -p genome && cd genome
# download Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
#      and Saccharomyces_cerevisiae.R64-1-1.114.gtf.gz  from ensembl.org
gunzip *.gz && cd ..

# 1b. verify FASTA and GTF agree on chromosome names — silent zero counts
#     downstream come from here and from nowhere else
comm -13 <(grep ">" genome/*.fa | cut -d" " -f1 | tr -d ">" | sort -u) \
         <(cut -f1 genome/*.gtf | grep -v "^#" | sort -u)

# 2. reads — accessions in data/samplesheet.csv, from ENA PRJEB5348
#    into data/raw/

# 3. QC, trim, index, align
bash src/fastqc.sh
bash src/trim.sh
bash src/index_creation_star.sh
bash src/alignment.sh

# 4. counts and differential expression
cd src && Rscript count_matrix.R && Rscript run_deseq2.R
```

`src/alignment.sh` clears each sample's output prefix before running, so it is
safe to re-run after a partial failure. Everything upstream of it is not yet
idempotent (see [Known limitations](#current-limitations)).



## Current limitations

1. **All shell scripts use paths relative to `src/`** and must be run from
   there (`cd src && bash alignment.sh`). They take no arguments, so input and
   output locations are fixed.
2. **Independent-filtering threshold mismatch.** `lfcShrink()` inherits the
   `results()` call made at the default `alpha = 0.1`, while significance is
   reported at 0.05. The optimal filtering threshold differs slightly between
   the two, so a small number of borderline padj values would change if the
   contrast were recomputed at `alpha = 0.05` throughout.
3. **Per-tile sequence quality fails in 11 of 12 samples.** This is a 2014
   HiSeq 2000 flow-cell artefact, not a library problem, and no reads were
   removed on account of it. Whether the failures correlate with lane has not
   been checked.
4. **No batch term in the design.** The model is `~ condition` only. The
   original study's samples span multiple lanes; lane was not tested as a
   covariate.

---


## References

- Gierliński M, *et al.* (2015). Statistical models for RNA-seq data derived from
  a two-condition 48-replicate experiment. *Bioinformatics* 31(22):3625–3630.
  [doi:10.1093/bioinformatics/btv425](https://doi.org/10.1093/bioinformatics/btv425)
- Schurch NJ, *et al.* (2016). How many biological replicates are needed in an
  RNA-seq experiment and which differential expression tool should you use?
  *RNA* 22(6):839–851. [doi:10.1261/rna.053959.115](https://doi.org/10.1261/rna.053959.115)
- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and
  dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550.
  [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
- Zhu A, Ibrahim JG, Love MI (2019). Heavy-tailed prior distributions for
  sequence count data (`apeglm`). *Bioinformatics* 35(12):2084–2092.
  [doi:10.1093/bioinformatics/bty895](https://doi.org/10.1093/bioinformatics/bty895)
- Dobin A, *et al.* (2013). STAR: ultrafast universal RNA-seq aligner.
  *Bioinformatics* 29(1):15–21.
  [doi:10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)

## License

MIT — see [LICENSE](LICENSE).
