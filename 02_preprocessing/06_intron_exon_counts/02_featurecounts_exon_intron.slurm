#!/bin/bash
# Counts exonic and intronic fragments for the embryo samples from the existing
# STAR alignments. Exon counts use every read pair overlapping a gene's exons;
# intron counts use only pairs lying wholly inside intron segments.
# NOTE: #SBATCH log paths are relative to the submission directory, which must
# contain a log/ folder (mkdir -p log before sbatch).
#SBATCH --job-name=fc_exon_intron
#SBATCH -o log/featurecounts_exon_intron_%j.out
#SBATCH -e log/featurecounts_exon_intron_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=08:00:00

# ---- Paths (set PROJECT_ROOT and ENV_PREFIX) --------------------------------
PROJECT_ROOT="/gpfs/home/jah14xyu/Postdoc/"
ENV_PREFIX="/gpfs/home/jah14xyu/micromamba/envs/fc_exon_intron"
SCRIPT_DIR="${PROJECT_ROOT}/ccap-rnaseq-organ-expression/02_preprocessing/06_intron_exon_counts"
BAM_DIR="${PROJECT_ROOT}/Bulk_RNA/results/results_short_read_genome/star_salmon"
BAM_SUFFIX="markdup.sorted.bam"
GTF="${PROJECT_ROOT}/Genomes/Short_Read/short_read_ccap21/GCF_000347755.3_genomic.fixed.gtf"
SAMPLE_LIST="${SCRIPT_DIR}/embryo_samples.txt"
OUT_DIR="${SCRIPT_DIR}/results/exon_intron_counts"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
MIN_INTRON_BP=100

# ---- Environment (created from environment.yml) ----------------------------
module purge
module load mamba
source activate "$ENV_PREFIX"
set -euo pipefail

# Env tools first on PATH; ignore R libraries outside the env.
export PATH="${ENV_PREFIX}/bin:${PATH}"
unset R_LIBS R_LIBS_USER

# ---- Upfront input checks --------------------------------------------------
command -v featureCounts >/dev/null || { echo "ERROR: featureCounts not found in $ENV_PREFIX (see environment.yml)" >&2; exit 1; }
Rscript -e 'library(rtracklayer); library(GenomicRanges); library(data.table)' >/dev/null 2>&1 \
    || { echo "ERROR: R packages missing in $ENV_PREFIX (see environment.yml)" >&2; exit 1; }
[[ -s "$GTF" ]]         || { echo "ERROR: GTF not found: $GTF" >&2; exit 1; }
[[ -s "$SAMPLE_LIST" ]] || { echo "ERROR: sample list not found: $SAMPLE_LIST" >&2; exit 1; }

BAMS=()
while read -r sample; do
    [[ -z "$sample" ]] && continue
    bam="${BAM_DIR}/${sample}.${BAM_SUFFIX}"
    [[ -s "$bam" ]] || { echo "ERROR: BAM not found: $bam (set BAM_SUFFIX if the STAR output is named differently)" >&2; exit 1; }
    BAMS+=("$bam")
done < "$SAMPLE_LIST"
echo "Samples: ${#BAMS[@]}"

mkdir -p "$OUT_DIR"

# ---- Step 1: exon and intron annotations -----------------------------------
Rscript "${SCRIPT_DIR}/01_build_exon_intron_features.R" "$GTF" "$OUT_DIR" "$MIN_INTRON_BP"

# ---- Step 2: exon counts ---------------------------------------------------
if [[ -s "${OUT_DIR}/exon_counts.txt" ]]; then
    echo "[OK] exon_counts.txt already exists - skipping"
else
    echo "[$(date)] Counting exonic fragments..."
    featureCounts -p --countReadPairs -B -C -s 0 -T "$THREADS" \
        -F SAF -a "${OUT_DIR}/exons.saf" \
        -o "${OUT_DIR}/exon_counts.tmp" "${BAMS[@]}"
    mv "${OUT_DIR}/exon_counts.tmp" "${OUT_DIR}/exon_counts.txt"
    mv "${OUT_DIR}/exon_counts.tmp.summary" "${OUT_DIR}/exon_counts.txt.summary"
fi

# ---- Step 3: intron counts -------------------------------------------------
if [[ -s "${OUT_DIR}/intron_counts.txt" ]]; then
    echo "[OK] intron_counts.txt already exists - skipping"
else
    echo "[$(date)] Counting intronic fragments..."
    featureCounts -p --countReadPairs -B -C -s 0 --fracOverlap 1 -T "$THREADS" \
        -F SAF -a "${OUT_DIR}/introns.saf" \
        -o "${OUT_DIR}/intron_counts.tmp" "${BAMS[@]}"
    mv "${OUT_DIR}/intron_counts.tmp" "${OUT_DIR}/intron_counts.txt"
    mv "${OUT_DIR}/intron_counts.tmp.summary" "${OUT_DIR}/intron_counts.txt.summary"
fi

echo "[OK] Outputs in ${OUT_DIR}:"
ls -la "${OUT_DIR}"/exon_counts.txt* "${OUT_DIR}"/intron_counts.txt* "${OUT_DIR}"/intron_lengths.tsv
