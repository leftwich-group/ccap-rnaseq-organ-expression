#!/bin/bash

set -euo pipefail

PROJECT_ROOT="/path/to/your/project"   # set this to your project root

# =============================================================================
# INPUT RBH FILES
# =============================================================================

SHORT_VS_LONG="${PROJECT_ROOT}/Bulk_RNA/results/making_gene_mapping_file/ccap21_vs_egii321/05_results/gene_mapping_rbh_filtered.tsv"
LONG_VS_DESC="${PROJECT_ROOT}/Bulk_RNA/results/making_gene_mapping_file/egii321_vs_deschepper/05_results/gene_mapping_rbh_filtered.tsv"

# =============================================================================
# OUTPUT
# =============================================================================

OUTDIR="${PROJECT_ROOT}/Bulk_RNA/results/making_gene_mapping_file/final_merged/merged_gene_map"
mkdir -p "${OUTDIR}"

# =============================================================================
# TEMP FILES
# =============================================================================

TMP1="${OUTDIR}/egii321_to_ccap21.tsv"
TMP2="${OUTDIR}/egii321_to_deschepper.tsv"

FINAL="${OUTDIR}/medfly_master_gene_map.tsv"

# =============================================================================
# EXTRACT egii321 ↔ ccap21
#
# short_vs_long structure:
#   short_gene_id = ccap21
#   long_gene_id  = egii321
#
# We reorder to:
#   egii321    ccap21
# =============================================================================

echo "[INFO] Processing short_vs_long..."

tail -n +2 "${SHORT_VS_LONG}" \
| cut -f1,7 \
| awk 'BEGIN{OFS="\t"}{print $2,$1}' \
| sort -k1,1 \
> "${TMP1}"

# =============================================================================
# EXTRACT egii321 ↔ deschepper
#
# long_vs_deschepper structure:
#   short_gene_id = egii321
#   long_gene_id  = deschepper
#
# Already ordered correctly:
#   egii321    deschepper
# =============================================================================

echo "[INFO] Processing long_vs_deschepper..."

tail -n +2 "${LONG_VS_DESC}" \
| cut -f1,7 \
| sort -k1,1 \
> "${TMP2}"

# =============================================================================
# JOIN TABLES
# =============================================================================

echo "[INFO] Joining tables..."

join \
    -t $'\t' \
    -a1 -a2 \
    -e NA \
    -o auto \
    "${TMP1}" \
    "${TMP2}" \
> "${OUTDIR}/joined.tmp"

# =============================================================================
# ADD HEADER
# =============================================================================

{
    echo -e "egii321_gene\tccap21_gene\tdeschepper_gene"
    cat "${OUTDIR}/joined.tmp"
} > "${FINAL}"

# =============================================================================
# CLEANUP
# =============================================================================

rm -f "${TMP1}" "${TMP2}" "${OUTDIR}/joined.tmp"

# =============================================================================
# SUMMARY
# =============================================================================

echo "[INFO] Done."
echo "[INFO] Output:"
echo "${FINAL}"

echo "[INFO] Number of merged genes:"
tail -n +2 "${FINAL}" | wc -l

