#!/bin/bash
# NOTE: #SBATCH directives are parsed literally by sbatch and cannot expand
# shell variables, so the log paths below are relative to the submission
# directory (expected to be $PROJECT_ROOT/Bulk_RNA, set below). Adjust to an
# absolute path if you submit this job from elsewhere.
#SBATCH --job-name=nfcore_rnaseq
#SBATCH -o log/RNAseq_nfcore_%j.out
#SBATCH -e log/RNAseq_nfcore_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@institution.ac.uk  # set your email here

set -euo pipefail

########################################
# 1) LOAD MODULES
########################################
module purge
module load nextflow
module load apptainer/1.4.0-rc.1
module load python/anaconda/2023.07/3.11

########################################
# 2) PATHS
########################################
PROJECT_ROOT="/path/to/your/project"   # set this to your project root
PROJECT_DIR="${PROJECT_ROOT}/Bulk_RNA"
RAW_DIR="${PROJECT_DIR}/raw_data/D:X204SC22083995-Z01-F001/01.RawData"
DATA_DIR="${PROJECT_DIR}/data"
LOG_DIR="${PROJECT_DIR}/log"
WORK_DIR="${PROJECT_DIR}/nextflow_work"
RESULTS_DIR="${PROJECT_DIR}/results/results_deschepper_genome"

SAMPLESHEET="${DATA_DIR}/samplesheet.csv"
CONFIG="${PROJECT_DIR}/data/rnaseq_slurm_deschepper_genome.config"

STRANDEDNESS="unstranded"

mkdir -p "${DATA_DIR}" "${LOG_DIR}" "${WORK_DIR}" "${RESULTS_DIR}"
mkdir -p "${PROJECT_DIR}/apptainer_cache" "${PROJECT_DIR}/tmp"

########################################
# 3) NEXTFLOW / APPTAINER CACHE LOCATIONS
########################################
export NXF_HOME="${PROJECT_DIR}/.nextflow"
export NXF_WORK="${WORK_DIR}"
export NXF_APPTAINER_CACHEDIR="${PROJECT_DIR}/apptainer_cache"
export APPTAINER_CACHEDIR="${PROJECT_DIR}/apptainer_cache"
export NXF_TEMP="${PROJECT_DIR}/tmp"
export TMPDIR="${PROJECT_DIR}/tmp"

########################################
# 4) BUILD SAMPLE SHEET
########################################
echo "sample,fastq_1,fastq_2,strandedness" > "${SAMPLESHEET}"

find "${RAW_DIR}" -mindepth 1 -maxdepth 1 -type d | sort | while read -r sample_dir; do
    sample=$(basename "${sample_dir}")

    r1="${sample_dir}/${sample}_1.fq.gz"
    r2="${sample_dir}/${sample}_2.fq.gz"

    if [[ ! -f "${r1}" ]]; then
        r1="${sample_dir}/${sample}_1.fastq.gz"
    fi
    if [[ ! -f "${r2}" ]]; then
        r2="${sample_dir}/${sample}_2.fastq.gz"
    fi

    if [[ ! -f "${r1}" ]]; then
        echo "ERROR: R1 not found for sample ${sample}" >&2
        exit 1
    fi
    if [[ ! -f "${r2}" ]]; then
        echo "ERROR: R2 not found for sample ${sample}" >&2
        exit 1
    fi

    r1_abs=$(readlink -f "${r1}")
    r2_abs=$(readlink -f "${r2}")

    echo "${sample},${r1_abs},${r2_abs},${STRANDEDNESS}" >> "${SAMPLESHEET}"
done

echo "Built samplesheet:"
cat "${SAMPLESHEET}"

########################################
# 5) RUN NF-CORE/RNASEQ
########################################
cd "${PROJECT_DIR}"

nextflow run nf-core/rnaseq \
    -r 3.23.0 \
    -profile apptainer \
    -c "${CONFIG}" \
    --input "${SAMPLESHEET}" \
    --fasta ${PROJECT_ROOT}/Genomes/Deschepper_2024/Ccapitata_DDBJ_100Ngaps.fa \
    --gff ${PROJECT_ROOT}/Genomes/Deschepper_2024/Ccapitata_DDBJ_100Ngaps.name_fixed.gff \
    --outdir ${PROJECT_ROOT}/Bulk_RNA/results/results_deschepper_genome \
    --email your.email@institution.ac.uk \
    --save_reference \
    -resume
