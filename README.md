# Ceratitis capitata bulk RNA-seq — organ-level expression atlas

Bulk RNA-seq preprocessing and downstream differential-expression analysis
of the Mediterranean fruit fly (medfly, *Ceratitis capitata*) across organs,
sexes, and developmental/reproductive states. This repository covers the
full path from raw reads to the biological results: read alignment and
quantification, cross-genome gene-ID reconciliation, GO annotation, and the
R analysis addressing four biological questions.

## Biological themes

The analysis is organised around four questions, plus a broader tissue atlas:

1. **Maternal-to-zygotic transition (MZT)** — transcriptional changes
   between unfertilised and fertilised eggs across early developmental
   timepoints.
2. **Sexual dimorphism** — sex-biased gene expression across germline
   (testes/ovary) and soma (abdomen) tissues.
3. **Post-mating transcriptional remodelling** — transcriptional response to
   mating in female reproductive and somatic tissues.
4. **Transcriptional determinants of phenotypic maleness** — expression
   signatures distinguishing male-specific tissues and developmental states.
5. **Tissue landscapes / global atlas** — a broader, cross-tissue view:
   PCA/UMAP of all tissues and sexes, tissue-specificity (tau index),
   constitutive vs. tissue-specific expression, and cross-tissue
   germline–soma correlation.

A companion analysis of Y-chromosome candidate genes and tephritid synteny,
which was originally interleaved with this analysis, has been split out
into a separate repository (`ccap-y-chromosome-synteny`) and is not
included here.

## Pipeline flow

```
Reference genome download
        │
        ▼
nf-core/rnaseq (read QC, trimming, STAR alignment, Salmon quantification)
        │
        ▼
Reciprocal-best-hit (RBH) gene-ID crosswalk across genome builds
        │
        ▼
GO annotation (eggNOG-mapper primary; InterProScan alternative)
        │
        ▼
R analysis: differential expression, GO enrichment, tissue-specificity,
            cross-tissue comparisons (03_analysis/Ccap_Bulk_Transcriptomics.Rmd)
```

## Reference genomes

Three independently assembled/annotated genome builds are used throughout,
which is why a gene-ID crosswalk step is needed before results from
different builds can be compared:

| Build | Accession | Description |
|---|---|---|
| Short-read | `GCF_000347755.3` (Ccap_2.1) | NCBI RefSeq, the most completely annotated reference; used as the primary genome for most analyses |
| Long-read | `GCA_905071925.1` (EGII-3.2.1) | GenBank long-read assembly |
| Deschepper 2024 | Zenodo record [13928607](https://zenodo.org/records/13928607) | Deschepper et al. 2024 assembly and annotation |

## Data availability

Raw FASTQ reads and per-sample BAM alignment files are **not included** in
this repository due to their size (BAM files run several GB per genome
build). To reproduce the count matrices from raw reads:

1. Download the reference genomes using the accessions above
   (`02_preprocessing/01_genome_download/`).
2. Run the nf-core/rnaseq pipeline against your own raw FASTQ files
   (`02_preprocessing/02_nfcore_rnaseq/`).

The small, final gene-level count and TPM matrices produced by this
pipeline (short-read and long-read genome builds) are included under
`02_preprocessing/05_results_summary/`.

A handful of large per-contrast R result objects (`results_all_nested.rds`,
tens of MB each, containing the full nested differential-expression result
list for a contrast) are likewise excluded from `04_results/` for size —
regenerate them by re-running the corresponding contrast section of
`03_analysis/Ccap_Bulk_Transcriptomics.Rmd`. All final summary tables,
figures, and TSV differential-expression outputs for every contrast are
included.

Raw sequencing reads are deposited at the European Nucleotide Archive (ENA)
under study accession [PRJEB124423](https://www.ebi.ac.uk/ena/browser/view/PRJEB124423).
BAM files and other large intermediates not covered above are available on
request from the repository maintainer.

## Repository layout

- **`02_preprocessing/`** — genome download, nf-core/rnaseq launch scripts and
  configs, the RBH gene-ID crosswalk pipeline, GO annotation guides, and
  summary count matrices. See `02_preprocessing/README.md`.
- **`03_analysis/`** — the R Markdown analysis (`Ccap_Bulk_Transcriptomics.Rmd`)
  covering all four biological questions plus the tissue atlas. See
  `03_analysis/README.md`.
- **`01_data/`** — small reference tables consumed by the R analysis: sample
  metadata, a biomaRt-style gene annotation object, and a digitised table
  of published sex-determining genes.
- **`04_results/`** — final output tables and figures from the R analysis,
  organised by theme (`01_MZT/` … `10_CrossTissue_correlation/`),
  plus `11_developmental_confirmation/` (QC figures) and
  `12_manuscript_plots/`.
- **`05_docs/`** — a narrative write-up of the preprocessing workflow, and a
  teaching walkthrough of the RNA-seq pipeline run manually on a single
  sample. See `05_docs/README.md`.

## Software

The RNA-seq quantification pipeline is
[nf-core/rnaseq](https://nf-co.re/rnaseq) v3.23.0, run via Nextflow with
Apptainer/Singularity containers under SLURM. Gene-ID crosswalking uses
DIAMOND blastp; GO annotation uses eggNOG-mapper (primary) and InterProScan 6
(alternative). The R analysis uses `edgeR`/`limma`-style quasi-likelihood
differential expression, `topGO` for GO enrichment, and standard
tidyverse/ggplot2 plotting. See `03_analysis/Ccap_Bulk_Transcriptomics.Rmd`
for the full package list.

## Running the pipeline yourself

Scripts under `02_preprocessing/` use a `PROJECT_ROOT` placeholder variable
in place of a real HPC path — set it to your own working directory before
running. `#SBATCH` directives (which cannot expand shell variables) use
relative log paths instead; adjust as needed for your submission setup.
`02_preprocessing/02_nfcore_rnaseq/samplesheet.csv` retains its original
raw-data delivery folder path (a Novogene sequencing order code) as
provenance; it will not resolve outside the cluster the pipeline was
originally run on.

A few `here()` calls in `03_analysis/Ccap_Bulk_Transcriptomics.Rmd` (the
eggNOG-mapper raw annotation input, the uncompressed reference GFF/FASTA used
by the HOMER motif chunk, and one legacy gene-map path) point at large or raw
files that are not included in this repository (see Data availability
above). Re-run those chunks after placing the corresponding files locally at
the referenced paths, or edit the paths to point at your own copies.
