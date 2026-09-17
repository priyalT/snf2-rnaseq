# snf2-rnaseq
## Summary

Bulk RNA-seq from raw reads: differential gene expression in *Saccharomyces cerevisiae* involving two conditions:
**Δsnf2** versus **wild type**.

Twelve samples were included (6 WT, 6 Δsnf2) from [ENA PRJEB5348](https://www.ebi.ac.uk/ena/browser/view/PRJEB5348)
taken end to end from: 

**raw FASTQ** → QC → adapter/quality trimming → splice-aware
alignment → gene-level counting → **differential expression**


## Result

The most significant change in the whole experiment is
**SNF2 itself**, the gene that was knocked out.

| Gene | Systematic ID | log2FC (shrunken) | baseMean | padj |
|---|---|---:|---:|---:|
| **SNF2** | YOR290C | **−5.35** | 120 | 1.5 × 10⁻⁵⁸ |

Summary of the full contrast (Δsnf2 vs WT, Wald test, Benjamini–Hochberg):

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
| `plots/PCA_plot.png` | VST-transformed samples; WT and Δsnf2 separate on PC1 |
| `plots/MA_plot.png` | Unshrunken log2 fold changes |
| `plots/MA_shrunkenLFC_plot.png` | The same contrast after `apeglm` shrinkage |
| `plots/ShrunkenLFC_comparison_plot.png` | Unshrunken vs shrunken LFC, with y = x |
| `plots/Volcano_plot.png` | padj < 0.05 & \|log2FC\| > 1 highlighted |



## Data

[Gierliński et al. 2015](https://doi.org/10.1093/bioinformatics/btv425) /
[Schurch et al. 2016](https://doi.org/10.1261/rna.053959.115) generated 48
biological replicates each of wild-type and Δsnf2 *S. cerevisiae*.

| | |
|---|---|
| Accession | ENA PRJEB5348 |
| Organism | *S. cerevisiae*, BY4741 background |
| Samples used | 12 of 96 (6 WT, 6 Δsnf2) — see `data/samplesheet.csv` |
| Reads | 51 bp, single-end, Illumina HiSeq 2000 (2014) |
| Reference | Ensembl R64-1-1, `dna.toplevel.fa` + release-114 GTF |

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

**Requirements:** `fastqc`, `multiqc`, `trim_galore` (+ `cutadapt`), `STAR ≥ 2.7`,
`R ≥ 4.3` with `DESeq2`, `apeglm`, `tidyverse`, `glue`, `readr`. Python side is
managed with [uv](https://docs.astral.sh/uv/) — `uv sync`.

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

1. **Strandedness was never verified.** Column 2 (unstranded) of
   `ReadsPerGene.out.tab` was used without any proper confirmation.
2. **`src/fastqc.sh` and `src/trim.sh` glob the current working directory** and
   take no arguments, so they only work if invoked from `data/raw/`. They are a
   record of what was run, not a runnable interface.
3. **`src/alignment.sh` hardcodes an absolute path** to the samplesheet.
4. **Independent-filtering threshold mismatch.** `lfcShrink()` inherits the
   `results()` call made at the default `alpha = 0.1`, while significance is
   reported at 0.05. The optimal filtering threshold differs slightly between
   the two, so a small number of borderline padj values would change if the
   contrast were recomputed at `alpha = 0.05` throughout.
7. **Per-tile sequence quality fails in 11 of 12 samples.** This is a 2014
   HiSeq 2000 flow-cell artefact, not a library problem, and no reads were
   removed on account of it. Whether the failures correlate with lane has not
   been checked.
8. **No batch term in the design.** The model is `~ condition` only. The
   original study's samples span multiple lanes; lane was not tested as a
   covariate.

---

## Future improvements

- [ ] Verify strandedness; re-run if column 2 was the wrong choice
- [ ] Make `fastqc.sh` / `trim.sh` / `alignment.sh` take paths as arguments
- [ ] Move package installation out of `run_deseq2.R` into an environment file
      (`environment.yml` / `renv.lock`)
- [ ] GO and pathway enrichment on the DE list — formally test the phosphate signal
- [ ] Lane as a covariate: `~ lane + condition`, compare
- [ ] **Replicate-depth analysis** — subsample *n* = 2, 3, 6, 12, 24, 48 per
      condition from the full 96-sample design, repeat the contrast, and plot how
      the recovered DE set and its fold changes stabilise with *n*. This is the
      question the dataset was built for and the natural extension of this repo.
- [ ] **Wrap the pipeline in Nextflow DSL2** — one process per step
      (`FASTQC`, `TRIMGALORE`, `STAR_GENOMEGENERATE`, `STAR_ALIGN`, `MULTIQC`),
      driven by `data/samplesheet.csv`, using nf-core modules where they exist,
      with containers per process and `-resume` for restartability. The shell
      scripts here map onto processes almost one-to-one, which is the point of
      having written them this way first.

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
