# Analysis

The downstream R analysis — quality control, quasi-likelihood differential
expression (`edgeR`), GO enrichment (`topGO`), and the tissue-specificity /
atlas analyses — split into one script per section under `scripts/`, driven by
a run-all pipeline. It was previously a single ~6,500-line R Markdown document,
now kept for reference at `_archive/Ccap_Bulk_Transcriptomics.Rmd`.

All paths resolve relative to the **repository root** with `here::here(...)`, so
every script must be run with the repository root as the working directory
(`run_pipeline.sh` and `run_all.R` `cd` there for you). `Bulk_RNA.Rproj` opens
an RStudio session here; set the working directory to the parent (repo root)
before sourcing scripts.

## Layout

```
03_analysis/
├── run_pipeline.sh        # run every stage in order (bash; per-stage logs)
├── run_all.R              # same, as an R wrapper
├── R/
│   ├── init.R             # bootstrap: sourced first by every script
│   ├── setup.R            # packages, palettes, themes, helper functions, constants
│   ├── paths.R            # canonical 04_results/ directory map (RES)
│   └── load_data.R        # count matrices, metadata, GO annotations, append_de_stats()
├── scripts/
│   ├── 00_init.R          # truncate 04_results/DE_stats_master.tsv for a fresh run
│   ├── 01_mzt.R                     -> 04_results/01_MZT/
│   ├── 02_sexual_dimorphism.R       -> 04_results/02_Sexual_dimorphism/
│   ├── 03_postmating.R             -> 04_results/03_PostMating/
│   ├── 04_maleness.R              -> 04_results/04_Maleness/       (reads 02)
│   ├── 05_tissue_landscapes.R     -> 04_results/05_Tissue_landscapes/
│   ├── 06_global_atlas.R          -> 04_results/06_Global_atlas/
│   ├── 07_tau_specificity.R       -> 04_results/07_Tau_specificity/
│   ├── 08_upset_sex_bias.R        -> 04_results/08_UpSet_sex_bias/  (reads 02)
│   ├── 09_homer_motifs.R          -> 04_results/09_HOMER_motifs/    (reads 02, 03; opt-in)
│   ├── 10_crosstissue_correlation.R -> 04_results/10_CrossTissue_correlation/
│   ├── 11_developmental_confirmation.R -> 04_results/11_developmental_confirmation/
│   ├── 12_manuscript_plots.R      -> 04_results/12_manuscript_plots/ (reads 08)
│   ├── 99_de_stats_master.R       -> 04_results/DE_stats_master_wide.tsv + overview figure
│   ├── 14_mzt_intron_exon.R       -> 04_results/14_MZT_intron_exon/  (needs 02_preprocessing/06_intron_exon_counts; run by hand)
│   ├── 15_mzt_provisioning.R      -> 04_results/15_MZT_provisioning/ (run by hand)
│   ├── 16_postmating_genesets.R   -> 04_results/16_PostMating_genesets/ (reads 03; run by hand)
│   ├── 17_supp_go_tissues.R       -> 04_results/17_Supp_GO_tissues/ (reads 05; run by hand)
│   ├── 18_housekeeping_tests.R    -> 04_results/18_Housekeeping_tests/ (reads 05; run by hand)
│   ├── 19_new_figures.R           -> 04_results/19_New_figures/ (assembles Figure 3 and Figures S14-S17 from 14-18; run by hand)
│   └── go_annotation_prep.R       -> 02_preprocessing/04_go_annotation/geneID2GO_*.rds (opt-in)
└── _archive/
    └── Ccap_Bulk_Transcriptomics.Rmd   # the original monolithic document
```

Script numbers match the `04_results/NN_*` output folders. The pipeline's
**run order** differs from the numbering (embryo QC first; `04`/`08`/`12`/`99`
run after their inputs) and is defined in `run_pipeline.sh`.

## Running

```bash
bash 03_analysis/run_pipeline.sh              # run / resume
env FORCE=1 bash 03_analysis/run_pipeline.sh  # re-run every stage from scratch
```

The `env FORCE=1` prefix works in bash, zsh and tcsh; the bare
`FORCE=1 bash ...` form is bash-only and errors in tcsh (`FORCE=1: Command not
found`). In tcsh you can also `setenv FORCE 1` first, then `unsetenv FORCE`.

or `Rscript 03_analysis/run_all.R` (same behaviour). Each stage's console
output is written to `03_analysis/logs/<stage>.log`; independent stages
continue if one fails and the run exits non-zero if any did.

### Resume / skip

A stage that finishes cleanly writes `03_analysis/logs/<stage>.ok`, and the
runner **skips any stage that already has one**. So after a failure part-way
through, just re-run `run_pipeline.sh` — completed stages are skipped and it
picks up where it stopped.

- Redo one stage: delete its marker, e.g. `rm 03_analysis/logs/04_maleness.ok`,
  then re-run the pipeline (or run the script directly).
- Redo everything: `env FORCE=1 bash 03_analysis/run_pipeline.sh`, or
  `rm 03_analysis/logs/*.ok`.

### GO enrichment cache

Every `topGO` run (the bulk of the pipeline's wall time) is cached to
`04_results/.go_cache/`, keyed by its exact inputs (gene set + genome universe
+ ontology). A stage that *does* re-run therefore reuses prior GO results
instead of recomputing them — a re-run after fixing a downstream bug takes
minutes, not hours. `FORCE=1` (or `GO_FORCE=1`) bypasses the cache; delete
`04_results/.go_cache/` to clear it (do this after rebuilding the GO
annotations with `go_annotation_prep.R`).

### A single section on its own

```bash
Rscript 03_analysis/scripts/02_sexual_dimorphism.R   # from the repo root
```

`append_de_stats()` replaces any existing row for the same
section/genome/contrast, so re-running one section does not double-count in
`DE_stats_master.tsv`. Run `00_init.R` first only for a completely fresh master
table.

## How the scripts share state

Each script is standalone: `R/init.R` reloads the count matrices, metadata and
GO annotations every time (a few seconds). Cross-section inputs are passed
**through disk** — `04_results/<section>/.../results_all_nested.rds` and the GO
`*_terms_merged.tsv` files — so `04` and `08` read `02`'s outputs, `12` reads
`08`'s, and so on. The `results_all_nested.rds` objects are tens of MB each and
are not tracked in this repo; regenerate them by running the section that
writes them.

## Inputs

- `01_data/metadata.csv` — sample sheet. Experiment-2 abdomen libraries
  (`RN34`–`RN36`) are contaminated and are filtered out in `R/load_data.R`, so
  they are excluded from every section.
- `02_preprocessing/05_results_summary/{short,long}_read_genome/salmon.merged.gene_counts.tsv`
- `02_preprocessing/04_go_annotation/geneID2GO_{short,long}.rds`
  (rebuild with `scripts/go_annotation_prep.R` if missing)

## HOMER motif enrichment

`scripts/09_homer_motifs.R` needs a local HOMER install plus a genome FASTA/GFF
(setup notes are in the archived Rmd). It is excluded from the run-all
pipeline: configure HOMER, set `RUN_HOMER_DISCOVERY <- TRUE` near the top of
the script, then run it by hand. Without that, the script only re-plots any
HOMER `knownResults.txt` files already present under `04_results/09_HOMER_motifs/`.
