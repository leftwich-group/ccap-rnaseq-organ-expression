# Preprocessing

Scripts and small reference outputs covering the path from raw sequencing
reads to gene-level count matrices and GO annotation, ready for the R
analysis in `../03_analysis/`.

- **`01_genome_download/`** — downloads the three reference genome builds
  (short-read Ccap_2.1, long-read EGII-3.2.1, Deschepper et al. 2024) from
  NCBI and Zenodo.
- **`02_nfcore_rnaseq/`** — sample-sheet construction and the nf-core/rnaseq
  v3.23.0 launch scripts/configs for each genome build. `samplesheet.csv`
  retains its original Novogene delivery-folder path for provenance; the
  path won't resolve outside the original HPC environment.
- **`03_gene_id_crosswalk_rbh/`** — reciprocal-best-hit (RBH) DIAMOND
  protein search between each pair of genome builds, merged into a single
  three-way gene-ID crosswalk table (`results/medfly_master_gene_map.tsv`).
- **`04_go_annotation/`** — guides for the two GO annotation pipelines used
  (eggNOG-mapper, primary; InterProScan, alternative), plus the small
  derived `geneID2GO_{short,long}.rds` lookup tables used by `topGO` in the
  R analysis.
- **`05_results_summary/`** — the final merged gene-level count and TPM
  matrices (short-read and long-read genome builds) produced by
  nf-core/rnaseq. Per-sample BAM files are not included (large and
  regenerable from raw reads); the Deschepper-build count matrices are not
  included pending a re-run of that pipeline branch. Both are available on
  request from the repository maintainer in the meantime.

- **`06_intron_exon_counts/`** - recounts the embryo samples from the existing
  STAR alignments as exonic and intronic fragments per gene, for
  `../03_analysis/scripts/14_mzt_intron_exon.R`. See its README.

Scripts use a `PROJECT_ROOT` placeholder variable in place of a real HPC
absolute path — set it near the top of each script before running.
`#SBATCH` log-path directives are relative to the submission directory
(SLURM does not expand shell variables in `#SBATCH` lines).

See `../05_docs/01_pipeline_log.md` for a narrative walkthrough of how this
pipeline was actually run, including a couple of annotation quirks (missing
`gene_id` attributes, scaffold-name mismatches) that needed manual fixes.
