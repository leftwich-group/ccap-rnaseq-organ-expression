#!/usr/bin/env bash
set -euo pipefail

# Make a working directory
mkdir -p ~/refs/ceratitis_capitata
cd ~/refs/ceratitis_capitata

# =============================================================================
# 1) SHORT-READ / WELL-ANNOTATED REFERENCE
#    NCBI RefSeq assembly: GCF_000347755.3 (Ccap_2.1)
# =============================================================================

mkdir -p short_read_ccap21
cd short_read_ccap21

# Genome FASTA
wget -O GCF_000347755.3_genomic.fna.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/347/755/GCF_000347755.3_Ccap_2.1/GCF_000347755.3_Ccap_2.1_genomic.fna.gz"

# Gene annotation (GFF3)
wget -O GCF_000347755.3_genomic.gff.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/347/755/GCF_000347755.3_Ccap_2.1/GCF_000347755.3_Ccap_2.1_genomic.gff.gz"

# Optional: protein and transcript FASTA
wget -O GCF_000347755.3_protein.faa.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/347/755/GCF_000347755.3_Ccap_2.1/GCF_000347755.3_Ccap_2.1_protein.faa.gz"

wget -O GCF_000347755.3_rna.fna.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/347/755/GCF_000347755.3_Ccap_2.1/GCF_000347755.3_Ccap_2.1_rna.fna.gz"

cd ..

# =============================================================================
# 2) LONG-READ REFERENCE
#    GenBank assembly: GCA_905071925.1 (EGII-3.2.1)
# =============================================================================

mkdir -p long_read_egii321
cd long_read_egii321

# Genome FASTA
wget -O GCA_905071925.1_genomic.fna.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/071/925/GCA_905071925.1_EGII-3.2.1/GCA_905071925.1_EGII-3.2.1_genomic.fna.gz"

# Assembly report
wget -O GCA_905071925.1_assembly_report.txt \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/071/925/GCA_905071925.1_EGII-3.2.1/GCA_905071925.1_EGII-3.2.1_assembly_report.txt"

# GTF file
wget -O GCA_905071925.1_genomic.gtf.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/071/925/GCA_905071925.1_EGII-3.2.1/GCA_905071925.1_EGII-3.2.1_genomic.gtf.gz"

cd ..

# =============================================================================
# 3) UNZIP
# =============================================================================

gunzip -f short_read_ccap21/*.gz
gunzip -f long_read_egii321/*.gz

echo "Download complete."

# =============================================================================
# 4) Download Deschepper genome assembly and annotation from Zenodo
# =============================================================================


mkdir -p Deschepper_2024
cd Deschepper_2024

wget -c https://zenodo.org/records/13928607/files/Ccapitata_DDBJ_100Ngaps.fa
wget -c https://zenodo.org/records/13928607/files/Ccapitata_DDBJ_100Ngaps.gff