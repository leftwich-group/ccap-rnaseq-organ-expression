#!/bin/bash
#SBATCH --job-name=medfly_rbh_template
#SBATCH -o log/rbh_%j.out
#SBATCH -e log/rbh_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --partition=hmem
#SBATCH --qos=hmem
#SBATCH --mem=64G
#SBATCH --time=10:00:00

set -euo pipefail

# NOTE: #SBATCH directives above cannot expand shell variables, so the log
# paths are relative to the submission directory. PROJECT_ROOT below is used
# for all other paths in this script.
PROJECT_ROOT="/path/to/your/project"   # set this to your project root

# =============================================================================
# USER CONFIG (EDIT ONLY THIS SECTION PER RUN)
# =============================================================================

# Comparison name 
COMPARISON="egii321_vs_deschepper"

# Assembly A
A_NAME="egii321"
A_FASTA="${PROJECT_ROOT}/Genomes/Long_Read/long_read_egii321/GCA_905071925.1_genomic.fna.gz"
A_GFF="${PROJECT_ROOT}/Genomes/Long_Read/long_read_egii321/GCA_905071925.1_EGII-3.2.1_genomic.gff.gz"

# Assembly B
B_NAME="deschepper"
B_FASTA="${PROJECT_ROOT}/Genomes/Deschepper_2024/Ccapitata_DDBJ_100Ngaps.fa"
B_GFF="${PROJECT_ROOT}/Genomes/Deschepper_2024/Ccapitata_DDBJ_100Ngaps.name_fixed.gff"

# =============================================================================
# OUTPUT SETTINGS
# =============================================================================

OUTDIR="${PROJECT_ROOT}/Bulk_RNA/results/making_gene_mapping_file/${COMPARISON}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Filtering thresholds
MIN_PID=50
MIN_QCOV=50
MIN_SCOV=50
MAX_EVALUE="1e-10"

mkdir -p "${OUTDIR}"/{01_inputs,02_parsed,03_proteins,04_diamond,05_results,logs}
cd "${OUTDIR}"

# =============================================================================
# SOFTWARE CHECKS
# =============================================================================

echo "[INFO] Checking required tools..."

command -v python3 >/dev/null 2>&1 || { echo "[ERROR] python3 not found"; exit 1; }
command -v gffread >/dev/null 2>&1 || { echo "[ERROR] gffread not found"; exit 1; }
command -v diamond >/dev/null 2>&1 || { echo "[ERROR] diamond not found"; exit 1; }

echo "[INFO] python3: $(python3 --version)"
echo "[INFO] gffread: $(gffread --version 2>/dev/null || true)"
echo "[INFO] diamond: $(diamond version 2>/dev/null || true)"

# =============================================================================
# 1. COPY / DECOMPRESS INPUTS
# =============================================================================

echo "[INFO] Preparing input files..."

# Assembly A
if [[ "${A_FASTA}" == *.gz ]]; then
    gunzip -c "${A_FASTA}" > 01_inputs/${A_NAME}.fa
else
    cp "${A_FASTA}" 01_inputs/${A_NAME}.fa
fi

if [[ "${A_GFF}" == *.gz ]]; then
    gunzip -c "${A_GFF}" > 01_inputs/${A_NAME}.gff3
else
    cp "${A_GFF}" 01_inputs/${A_NAME}.gff3
fi

# Assembly B
if [[ "${B_FASTA}" == *.gz ]]; then
    gunzip -c "${B_FASTA}" > 01_inputs/${B_NAME}.fa
else
    cp "${B_FASTA}" 01_inputs/${B_NAME}.fa
fi

if [[ "${B_GFF}" == *.gz ]]; then
    gunzip -c "${B_GFF}" > 01_inputs/${B_NAME}.gff3
else
    cp "${B_GFF}" 01_inputs/${B_NAME}.gff3
fi

# =============================================================================
# 2. HELPER PYTHON SCRIPTS
# =============================================================================

cat > 02_parsed/parse_gff3_tx2gene.py <<'PY'
#!/usr/bin/env python3
import sys
import csv

gff_file = sys.argv[1]
out_file = sys.argv[2]

def parse_attr(attr_str):
    d = {}
    for item in attr_str.strip().split(';'):
        if not item:
            continue
        if '=' in item:
            k, v = item.split('=', 1)
            d[k] = v
        elif ' ' in item:
            # loose fallback for GTF-like attributes
            parts = item.strip().split(' ', 1)
            if len(parts) == 2:
                d[parts[0]] = parts[1].strip().strip('"')
    return d

gene_name_by_id = {}
gene_biotype_by_id = {}
tx_rows = []

tx_types = {
    "mRNA", "transcript", "RNA", "lnc_RNA", "ncRNA",
    "primary_transcript", "pseudogenic_transcript"
}

with open(gff_file) as fh:
    for line in fh:
        if not line.strip() or line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 9:
            continue

        seqid, source, feature_type, start, end, score, strand, phase, attrs = cols
        a = parse_attr(attrs)

        if feature_type == "gene":
            gene_id = a.get("ID") or a.get("gene_id")
            gene_name = a.get("Name") or a.get("gene") or a.get("gene_name") or ""
            biotype = (
                a.get("gene_biotype") or a.get("gene_type") or a.get("gbkey") or ""
            )
            if gene_id:
                gene_name_by_id[gene_id] = gene_name
                gene_biotype_by_id[gene_id] = biotype

        elif feature_type in tx_types:
            tx_id = a.get("ID") or a.get("transcript_id")
            parent = a.get("Parent") or a.get("gene_id")
            tx_name = a.get("Name") or a.get("transcript_name") or ""

            if tx_id and parent:
                gene_id = parent.split(',')[0]
                tx_rows.append((tx_id, gene_id, tx_name))

with open(out_file, 'w', newline='') as out:
    w = csv.writer(out, delimiter='\t')
    w.writerow(["tx_id", "gene_id", "gene_name", "tx_name", "gene_biotype"])
    for tx_id, gene_id, tx_name in tx_rows:
        w.writerow([
            tx_id,
            gene_id,
            gene_name_by_id.get(gene_id, ""),
            tx_name,
            gene_biotype_by_id.get(gene_id, "")
        ])
PY
chmod +x 02_parsed/parse_gff3_tx2gene.py

cat > 03_proteins/select_${B_NAME}est_isoform.py <<'PY'
#!/usr/bin/env python3
import sys
import csv

protein_fa = sys.argv[1]
tx2gene_tsv = sys.argv[2]
out_fa = sys.argv[3]
out_map = sys.argv[4]

# Protein symbols tolerated after cleaning.
# We strip terminal '*' before length comparison / output.
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")

tx2gene = {}
with open(tx2gene_tsv) as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        tx2gene[row["tx_id"]] = row

def read_fasta(path):
    name = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, ''.join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            yield name, ''.join(seq)

def clean_protein(seq):
    seq = seq.upper().rstrip('*')
    seq = seq.replace('.', 'X')
    seq = ''.join(aa if aa in VALID_AA else 'X' for aa in seq)
    return seq

best = {}  # gene_id -> record

for tx_id, seq in read_fasta(protein_fa):
    if tx_id not in tx2gene:
        continue

    seq = clean_protein(seq)

    if not seq:
        continue

    row = tx2gene[tx_id]
    gene_id = row["gene_id"]
    gene_name = row["gene_name"]
    tx_name = row["tx_name"]
    biotype = row["gene_biotype"]
    aa_len = len(seq)

    rec = {
        "gene_id": gene_id,
        "gene_name": gene_name,
        "tx_id": tx_id,
        "tx_name": tx_name,
        "gene_biotype": biotype,
        "aa_len": aa_len,
        "seq": seq
    }

    if gene_id not in best:
        best[gene_id] = rec
    else:
        old = best[gene_id]
        if aa_len > old["aa_len"]:
            best[gene_id] = rec
        elif aa_len == old["aa_len"] and tx_id < old["tx_id"]:
            best[gene_id] = rec

with open(out_fa, "w") as fa, open(out_map, "w", newline="") as tab:
    w = csv.writer(tab, delimiter='\t')
    w.writerow(["gene_id", "gene_name", "tx_id", "tx_name", "gene_biotype", "aa_len"])
    for gene_id in sorted(best):
        rec = best[gene_id]

        # IMPORTANT: FASTA header is gene_id only, so DIAMOND IDs match map keys.
        fa.write(f'>{rec["gene_id"]}\n')

        seq = rec["seq"]
        for i in range(0, len(seq), 80):
            fa.write(seq[i:i+80] + '\n')

        w.writerow([
            rec["gene_id"], rec["gene_name"], rec["tx_id"],
            rec["tx_name"], rec["gene_biotype"], rec["aa_len"]
        ])
PY
chmod +x 03_proteins/select_${B_NAME}est_isoform.py

cat > 05_results/rbh_from_diamond.py <<'PY'
#!/usr/bin/env python3
import sys
import csv

ab_file = sys.argv[1]
ba_file = sys.argv[2]
a_map_file = sys.argv[3]
b_map_file = sys.argv[4]
raw_out = sys.argv[5]
filt_out = sys.argv[6]
min_pid = float(sys.argv[7])
min_qcov = float(sys.argv[8])
min_scov = float(sys.argv[9])
max_evalue = float(sys.argv[10])

def load_map(path):
    d = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            d[row["gene_id"]] = row
    return d

def load_best_hits(path):
    best = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            qseqid, sseqid, bitscore, evalue, pident, qcovhsp, scovhsp = line.rstrip('\n').split('\t')
            rec = {
                "qseqid": qseqid,
                "sseqid": sseqid,
                "bitscore": float(bitscore),
                "evalue": float(evalue),
                "pident": float(pident),
                "qcovhsp": float(qcovhsp),
                "scovhsp": float(scovhsp),
            }
            if qseqid not in best:
                best[qseqid] = rec
            else:
                old = best[qseqid]
                if (
                    rec["bitscore"] > old["bitscore"] or
                    (rec["bitscore"] == old["bitscore"] and rec["evalue"] < old["evalue"]) or
                    (rec["bitscore"] == old["bitscore"] and rec["evalue"] == old["evalue"] and rec["pident"] > old["pident"])
                ):
                    best[qseqid] = rec
    return best

a_map = load_map(a_map_file)
b_map = load_map(b_map_file)

ab = load_best_hits(ab_file)
ba = load_best_hits(ba_file)

headers = [
    "${A_NAME}_gene_id", "${A_NAME}_gene_name", "${A_NAME}_tx_id", "${A_NAME}_tx_name", "${A_NAME}_gene_biotype", "${A_NAME}_aa_len",
    "${B_NAME}_gene_id", "${B_NAME}_gene_name", "${B_NAME}_tx_id", "${B_NAME}_tx_name", "${B_NAME}_gene_biotype", "${B_NAME}_aa_len",
    "pident", "qcovhsp", "scovhsp", "bitscore", "evalue", "method"
]

raw_rows = []
for a_gene, rec in ab.items():
    b_gene = rec["sseqid"]
    if b_gene in ba and ba[b_gene]["sseqid"] == a_gene:
        if a_gene not in a_map or b_gene not in b_map:
            continue
        a = a_map[a_gene]
        b = b_map[b_gene]
        raw_rows.append({
            "${A_NAME}_gene_id": a["gene_id"],
            "${A_NAME}_gene_name": a["gene_name"],
            "${A_NAME}_tx_id": a["tx_id"],
            "${A_NAME}_tx_name": a["tx_name"],
            "${A_NAME}_gene_biotype": a["gene_biotype"],
            "${A_NAME}_aa_len": a["aa_len"],
            "${B_NAME}_gene_id": b["gene_id"],
            "${B_NAME}_gene_name": b["gene_name"],
            "${B_NAME}_tx_id": b["tx_id"],
            "${B_NAME}_tx_name": b["tx_name"],
            "${B_NAME}_gene_biotype": b["gene_biotype"],
            "${B_NAME}_aa_len": b["aa_len"],
            "pident": f'{rec["pident"]:.2f}',
            "qcovhsp": f'{rec["qcovhsp"]:.2f}',
            "scovhsp": f'{rec["scovhsp"]:.2f}',
            "bitscore": f'{rec["bitscore"]:.2f}',
            "evalue": f'{rec["evalue"]:.3e}',
            "method": "diamond_blastp_rbh_${B_NAME}est_isoform"
        })

raw_rows = sorted(raw_rows, key=lambda x: (x["${A_NAME}_gene_id"], x["${B_NAME}_gene_id"]))

with open(raw_out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=headers, delimiter='\t')
    w.writeheader()
    for row in raw_rows:
        w.writerow(row)

with open(filt_out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=headers, delimiter='\t')
    w.writeheader()
    for row in raw_rows:
        if (
            float(row["pident"]) >= min_pid and
            float(row["qcovhsp"]) >= min_qcov and
            float(row["scovhsp"]) >= min_scov and
            float(row["evalue"]) <= max_evalue
        ):
            w.writerow(row)
PY
chmod +x 05_results/rbh_from_diamond.py

# =============================================================================
# 3. PARSE GFF3 TO TX->GENE TABLES
# =============================================================================

echo "[INFO] Parsing GFF3 files..."

python3 02_parsed/parse_gff3_tx2gene.py 01_inputs/${A_NAME}.gff3 02_parsed/${A_NAME}.tx2gene.tsv
python3 02_parsed/parse_gff3_tx2gene.py 01_inputs/${B_NAME}.gff3  02_parsed/${B_NAME}.tx2gene.tsv

# =============================================================================
# 4. EXTRACT PROTEINS
# =============================================================================

echo "[INFO] Extracting proteins with gffread..."

gffread 01_inputs/${A_NAME}.gff3 -g 01_inputs/${A_NAME}.fa -y 03_proteins/${A_NAME}.all_proteins.fa
gffread 01_inputs/${B_NAME}.gff3  -g 01_inputs/${B_NAME}.fa  -y 03_proteins/${B_NAME}.all_proteins.fa

# =============================================================================
# 5. KEEP ${B_NAME}EST ISOFORM PER GENE
# =============================================================================

echo "[INFO] Selecting ${B_NAME}est protein isoform per gene..."

python3 03_proteins/select_${B_NAME}est_isoform.py \
    03_proteins/${A_NAME}.all_proteins.fa \
    02_parsed/${A_NAME}.tx2gene.tsv \
    03_proteins/${A_NAME}.${B_NAME}est_gene_proteins.fa \
    03_proteins/${A_NAME}.${B_NAME}est_gene_proteins.tsv

python3 03_proteins/select_${B_NAME}est_isoform.py \
    03_proteins/${B_NAME}.all_proteins.fa \
    02_parsed/${B_NAME}.tx2gene.tsv \
    03_proteins/${B_NAME}.${B_NAME}est_gene_proteins.fa \
    03_proteins/${B_NAME}.${B_NAME}est_gene_proteins.tsv

# =============================================================================
# 5b. VALIDATE PROTEIN FASTA
# =============================================================================

echo "[INFO] Validating cleaned protein FASTA files..."

echo "[INFO] Checking for remaining dots in ${A_NAME} proteins..."
if grep -n '^[^>].*\.' 03_proteins/${A_NAME}.${B_NAME}est_gene_proteins.fa | head; then
    echo "[ERROR] Dots remain in ${A_NAME} protein FASTA."
    exit 1
else
    echo "[INFO] No dots found in ${A_NAME} protein FASTA."
fi

echo "[INFO] Checking for remaining dots in ${B_NAME} proteins..."
if grep -n '^[^>].*\.' 03_proteins/${B_NAME}.${B_NAME}est_gene_proteins.fa | head; then
    echo "[ERROR] Dots remain in ${B_NAME} protein FASTA."
    exit 1
else
    echo "[INFO] No dots found in ${B_NAME} protein FASTA."
fi

echo "[INFO] Checking for invalid amino-acid characters in ${A_NAME} proteins..."
awk '
  /^>/ {next}
  /[^ACDEFGHIKLMNPQRSTVWYBXZJUO]/ {
    print "[ERROR] Invalid ${A_NAME} protein line:", NR, $0
    exit 1
  }
' 03_proteins/${A_NAME}.${B_NAME}est_gene_proteins.fa

echo "[INFO] Checking for invalid amino-acid characters in ${B_NAME} proteins..."
awk '
  /^>/ {next}
  /[^ACDEFGHIKLMNPQRSTVWYBXZJUO]/ {
    print "[ERROR] Invalid ${B_NAME} protein line:", NR, $0
    exit 1
  }
' 03_proteins/${B_NAME}.${B_NAME}est_gene_proteins.fa

echo "[INFO] Protein FASTA validation passed."

# =============================================================================
# 6. BUILD DIAMOND DATABASES
# =============================================================================

echo "[INFO] Building DIAMOND databases..."

diamond makedb \
    --in 03_proteins/${A_NAME}.${B_NAME}est_gene_proteins.fa \
    --db 04_diamond/${A_NAME}.${B_NAME}est_gene_proteins

diamond makedb \
    --in 03_proteins/${B_NAME}.${B_NAME}est_gene_proteins.fa \
    --db 04_diamond/${B_NAME}.${B_NAME}est_gene_proteins

# =============================================================================
# 7. DIAMOND BLASTP BOTH WAYS
# =============================================================================

echo "[INFO] Running DIAMOND blastp..."

diamond blastp \
    --query 03_proteins/${A_NAME}.${B_NAME}est_gene_proteins.fa \
    --db 04_diamond/${B_NAME}.${B_NAME}est_gene_proteins.dmnd \
    --out 04_diamond/${A_NAME}_vs_${B_NAME}.tsv \
    --outfmt 6 qseqid sseqid bitscore evalue pident qcovhsp scovhsp \
    --max-target-seqs 5 \
    --evalue 1e-5 \
    --sensitive \
    --threads "${THREADS}"

diamond blastp \
    --query 03_proteins/${B_NAME}.${B_NAME}est_gene_proteins.fa \
    --db 04_diamond/${A_NAME}.${B_NAME}est_gene_proteins.dmnd \
    --out 04_diamond/${B_NAME}_vs_${A_NAME}.tsv \
    --outfmt 6 qseqid sseqid bitscore evalue pident qcovhsp scovhsp \
    --max-target-seqs 5 \
    --evalue 1e-5 \
    --sensitive \
    --threads "${THREADS}"

# =============================================================================
# 8. BUILD RBH GENE MAPPING TABLE
# =============================================================================

echo "[INFO] Building RBH gene mapping tables..."

python3 05_results/rbh_from_diamond.py \
    04_diamond/${A_NAME}_vs_${B_NAME}.tsv \
    04_diamond/${B_NAME}_vs_${A_NAME}.tsv \
    03_proteins/${A_NAME}.${B_NAME}est_gene_proteins.tsv \
    03_proteins/${B_NAME}.${B_NAME}est_gene_proteins.tsv \
    05_results/gene_mapping_rbh_raw.tsv \
    05_results/gene_mapping_rbh_filtered.tsv \
    "${MIN_PID}" "${MIN_QCOV}" "${MIN_SCOV}" "${MAX_EVALUE}"

# =============================================================================
# 9. SIMPLE SUMMARIES
# =============================================================================

echo "[INFO] Writing summaries..."

{
    echo -e "file\tcount"
    echo -e "${A_NAME}_tx2gene\t$(tail -n +2 02_parsed/${A_NAME}.tx2gene.tsv | wc -l)"
    echo -e "${B_NAME}_tx2gene\t$(tail -n +2 02_parsed/${B_NAME}.tx2gene.tsv | wc -l)"
    echo -e "${A_NAME}_${B_NAME}est_genes\t$(tail -n +2 03_proteins/${A_NAME}.${B_NAME}est_gene_proteins.tsv | wc -l)"
    echo -e "${B_NAME}_${B_NAME}est_genes\t$(tail -n +2 03_proteins/${B_NAME}.${B_NAME}est_gene_proteins.tsv | wc -l)"
    echo -e "rbh_raw\t$(tail -n +2 05_results/gene_mapping_rbh_raw.tsv | wc -l)"
    echo -e "rbh_filtered\t$(tail -n +2 05_results/gene_mapping_rbh_filtered.tsv | wc -l)"
} > 05_results/mapping_summary.tsv

echo "[INFO] Done."
echo "[INFO] Main output: ${OUTDIR}/05_results/gene_mapping_rbh_filtered.tsv"

# ====== END OF medfly_rbh_map.slurm ======

