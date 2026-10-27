# Exon and intron counts for the embryo time course

Recounts the embryo (egg) samples from the existing STAR alignments as exonic
and intronic fragments per gene, for the analysis in
`../../03_analysis/scripts/14_mzt_intron_exon.R`. No realignment is needed.

## Contents

- `01_build_exon_intron_features.R` - builds `exons.saf`, `introns.saf` and
  `intron_lengths.tsv` from the short-read (Ccap_2.1) GTF. Exons are the union
  across transcripts. Introns are the gene span minus the exons of every gene
  on either strand, because the libraries are unstranded.
- `02_featurecounts_exon_intron.slurm` - runs step 1 and two `featureCounts`
  jobs (exons, then introns) on the 27 embryo BAMs listed in
  `embryo_samples.txt`. Intron counts use `--fracOverlap 1`, so only fragments
  lying wholly inside an intron are counted.
- `results/` - the small tables that the R analysis reads (see below).

## Environment

`featureCounts` and the R packages come from a mamba environment. Create it
once on the HPC from `environment.yml`:

```bash
module load mamba
mamba env create -f environment.yml
```

Then set `ENV_PREFIX` in the slurm script to the new environment's path (for
example `/gpfs/home/<user>/micromamba/envs/fc_exon_intron`). Print the path
with `mamba env list`.

## Run

1. Set `PROJECT_ROOT` and `ENV_PREFIX` at the top of the slurm script. `BAM_SUFFIX` is the
   suffix of the nf-core STAR output (`markdup.sorted.bam` by default).
2. Submit from a directory that contains a `log/` folder:

```bash
mkdir -p log
sbatch 02_featurecounts_exon_intron.slurm
```

### Running step 1 in RStudio or Positron

Open `01_build_exon_intron_features.R`, check the three settings in the
interactive branch (annotation file, output folder, minimum intron length),
and run or source the file. The defaults read the repository's Ccap_2.1 GFF3
and write to `results/`. This is only for inspecting the annotation. The slurm
job rebuilds the files from the GTF that was used for alignment.

Only step 1 can run locally. The `featureCounts` steps need the STAR BAM files,
which are on the HPC, so run `02_featurecounts_exon_intron.slurm` there.

Steps skip if their outputs already exist (delete a file to redo that step).

3. Copy these files from `OUT_DIR` into `results/`:
   `exon_counts.txt`, `intron_counts.txt`, `intron_lengths.tsv`, and the two
   `.summary` files.
