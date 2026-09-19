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

The largest fold change in the whole experiment is **SNF2 itself**, the gene
that was knocked out. Of 1,372 significant genes, none moves further than the
one that was deleted.

| Gene | Systematic ID | log2FC (shrunken) | lfcSE | baseMean | padj |
|---|---|---:|---:|---:|---:|
| **SNF2** | YOR290C | **−6.92** | 0.42 | 120 | 1.6 × 10⁻⁵⁸ |

Summary of the full contrast (snf2Δ vs WT, Wald test, Benjamini–Hochberg):

| | |
|---|---:|
| Gene features in the annotation (Ensembl 114 on R64-1-1; 6,600 protein-coding) | 7,127 |
| Genes passing the count filter (`rowSums(counts) >= 10`) | 6,115 |
| Genes with a usable padj | 5,844 |
| Significant at padj < 0.05 | 1,372 |
| Significant **and** \|log2FC\| > 1 | 396 (119 up, 277 down) |

A quarter of the tested transcriptome moves at padj < 0.05. That is expected because
Snf2 is the catalytic ATPase of the SWI/SNF chromatin-remodelling complex, so
deleting it perturbs transcription globally rather than in one pathway. 

**Observed pattern:** the second and third largest effects after SNF2 are
`PHO12` and `SPL2`, and `PHO84`, `VTC3` and `VTC4` are close behind — all
phosphate-responsive genes, all down.
SWI/SNF is known to be required for chromatin remodelling at PHO promoters, so
this is a coherent signal rather than noise, and it was not something the
analysis was set up to look for.

### Statistical choices

**`alpha = 0.05` is passed explicitly to `results()`, and that result is passed
into `lfcShrink()`.** This matters more than it looks. DESeq2 does *independent
filtering* before BH correction — it discards genes whose mean normalised count
is too low to ever reach significance, which reduces the number of tests and so
reduces the multiple-testing penalty on the genes that remain. The threshold is
not fixed: DESeq2 chooses the quantile that maximises the number of genes
significant **at the alpha it was given**. Measured on this dataset:

| alpha | quantile chosen | baseMean threshold |
|---|---:|---:|
| 0.10 (the default) | 0.00 % | 0.64 |
| 0.05 (used here) | 3.88 % | 2.12 |

Running at the default 0.10 and then reporting at 0.05 changed the significance
call for 20 genes — 15 gained significance, 5 lost it. The 15 gains are the
point of the filter: removing ~240 hopeless genes lightens BH's penalty enough
for borderline genes to cross. The 5 losses all had `baseMean` below 2.12 and
nominal p-values around 0.005 on one or two reads per sample, which is exactly
what the filter exists to discard.

**The 271 genes with `padj = NA` come from two different mechanisms**, which are
distinguishable by whether `pvalue` is also `NA`:

- **238 independent filtering** — `pvalue` present, `padj` `NA`, all with
  `baseMean` ≤ 2.117, the threshold above.
- **33 Cook's-distance outliers** — both `pvalue` and `padj` `NA`. At least one
  sample carries a count extreme enough to dominate the fit. With 6 replicates
  per condition DESeq2 is below its `minReplicatesForReplace = 7` threshold, so
  it declines to test these rather than replacing the outlying count.

**Shrunken fold changes depend on the gene set they were fitted alongside.**
`apeglm` estimates its prior empirically, from the spread of fold changes across
every gene in the object. An earlier run of this analysis used a count matrix
with 5,401 genes and reported SNF2 at −5.35 (lfcSE 0.23); the current 6,115-gene
matrix gives −6.92 (lfcSE 0.42) for identical counts — `baseMean` 120.3 and
p = 1.4 × 10⁻⁶¹ in both. The extra low-count genes widened the observed fold-change
distribution, so the fitted prior is heavier-tailed and large effects are shrunk
less. A shrunken log2 fold change is therefore not a property of a gene on its
own, and the gene set is part of the result.

### Figures

| | |
|---|---|
| `plots/PCA_plot.png` | VST-transformed samples; WT and snf2Δ separate on PC1 |
| `plots/MA_plot.png` | Unshrunken log2 fold changes |
| `plots/MA_shrunkenLFC_plot.png` | The same contrast after `apeglm` shrinkage |
| `plots/compare_shrunkenLFC_plot.png` | Unshrunken vs shrunken LFC, with y = x |
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
2. **Per-tile sequence quality fails in 11 of 12 samples.** This is a 2014
   HiSeq 2000 flow-cell artefact, not a library problem, and no reads were
   removed on account of it. Whether the failures correlate with lane has not
   been checked.
3. **No batch term in the design.** The model is `~ condition` only. The
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
