# InterProScan gene-to-GO annotation guide

This is an alternative/cross-check GO annotation source (the primary source
used downstream is eggNOG-mapper — see `eggNOG-mapper_guide.md`). It
documents how a gene-to-GO mapping table was built from InterProScan 6 run
via Nextflow/Apptainer on the HPC cluster, for both the short-read and
long-read proteomes.

Throughout, `$PROJECT_ROOT` refers to your own working directory on the
cluster (e.g. `/path/to/your/project`).

## 1. Load modules

There is no direct HPC module for InterProScan; it is run via Nextflow and
Apptainer (Singularity-compatible).

```bash
module load nextflow/25.04.6
module load apptainer/1.4.0-rc.1

mkdir -p ${PROJECT_ROOT}/GO_Terms/interproscan
```

## 2. Clean the input protein FASTA

InterProScan is strict about its amino-acid alphabet, so the protein FASTA
needs to be cleaned first: whitespace stripped, stop codons (`*`) and gap
characters removed, and rare/ambiguous residues (`U`, `O`, `J`, `Z`) or
anything else outside the standard 20 amino acids + `B`/`X` mapped to `X`.

`clean_protein_fasta.sh` (a small Python-in-bash wrapper) does this:

```bash
#!/bin/bash
set -euo pipefail

IN_FASTA="$1"
OUT_FASTA="$2"
BAD_FASTA="$3"

python3 - <<EOF
import gzip
import re

infile = "${IN_FASTA}"
outfile = "${OUT_FASTA}"
badfile = "${BAD_FASTA}"

# Conservative amino-acid alphabet for InterProScan protein input
valid = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")

def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")

def flush_record(header, seq, out, bad):
    if header is None:
        return

    raw_seq = "".join(seq).upper()
    raw_seq = re.sub(r"\s+", "", raw_seq)

    cleaned = raw_seq
    cleaned = cleaned.replace("*", "").replace(".", "")
    cleaned = cleaned.replace("-", "").replace("_", "")

    # Convert rare/ambiguous residues some InterProScan member tools reject
    for ch in "UOJZ":
        cleaned = cleaned.replace(ch, "X")

    final = []
    bad_chars = set()
    for aa in cleaned:
        if aa in set("ACDEFGHIKLMNPQRSTVWYBX"):
            final.append(aa)
        else:
            final.append("X")
            bad_chars.add(aa)
    final = "".join(final)

    if len(final) == 0:
        bad.write(f"{header}\tEMPTY_AFTER_CLEANING\n")
        return

    if bad_chars or raw_seq != final:
        bad.write(f"{header}\tcleaned\tbad_chars={','.join(sorted(bad_chars))}\n")

    out.write(header + "\n")
    for i in range(0, len(final), 60):
        out.write(final[i:i+60] + "\n")

with open_maybe_gz(infile) as fh, open(outfile, "w") as out, open(badfile, "w") as bad:
    header = None
    seq = []
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush_record(header, seq, out, bad)
            header = line
            seq = []
        else:
            seq.append(line)
    flush_record(header, seq, out, bad)
EOF
```

Run it for both proteomes:

```bash
chmod +x clean_protein_fasta.sh

./clean_protein_fasta.sh \
  ${PROJECT_ROOT}/Genomes/Short_Read/short_read_ccap21/GCF_000347755.3_protein.faa.gz \
  ${PROJECT_ROOT}/GO_Terms/interproscan/data/proteomes/short_clean.faa \
  ${PROJECT_ROOT}/GO_Terms/interproscan/data/proteomes/short_cleaning_log.tsv

./clean_protein_fasta.sh \
  ${PROJECT_ROOT}/Genomes/Long_Read/long_read_egii321/GCA_905071925.1_EGII-3.2.1_protein.faa.gz \
  ${PROJECT_ROOT}/GO_Terms/interproscan/data/proteomes/long_clean.faa \
  ${PROJECT_ROOT}/GO_Terms/interproscan/data/proteomes/long_cleaning_log.tsv
```

Sanity check — this should print nothing:

```bash
awk '!/^>/' ${PROJECT_ROOT}/GO_Terms/interproscan/data/proteomes/short_clean.faa \
  | grep -n '[^ACDEFGHIKLMNPQRSTVWYBX]' | head

awk '!/^>/' ${PROJECT_ROOT}/GO_Terms/interproscan/data/proteomes/long_clean.faa \
  | grep -n '[^ACDEFGHIKLMNPQRSTVWYBX]' | head
```

## 3. Nextflow config

```groovy
process {
  executor = 'slurm'
  scratch = true
}

apptainer {
  enabled = true
  autoMounts = true
  cacheDir = "${PROJECT_ROOT}/.apptainer_cache"
}

singularity { enabled = false }
docker     { enabled = false }

workDir = "${PROJECT_ROOT}/GO_terms/interproscan"
```

## 4. Test run (50 sequences)

Before committing to a full run, test on a small subset:

```bash
LONG_FASTA="${PROJECT_ROOT}/Genomes/Long_Read/long_read_egii321/GCA_905071925.1_EGII-3.2.1_protein.faa.gz"

awk '/^>/{n++} n<=50 {print}' ${LONG_FASTA} > ${PROJECT_ROOT}/GO_Terms/interproscan/data/proteomes/test50.fa

module load nextflow/25.04.6
module load apptainer/1.4.0-rc.1

nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \
  -profile slurm,apptainer \
  --input data/proteomes/test50.fa \
  --datadir data/interpro_data \
  --outdir results/test50 \
  --formats tsv \
  --goterms \
  -resume
```

## 5. Full run (both proteomes)

`run_interproscan_long_short.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=interproscan_medfly_both
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=log/interpro_%x_%j.out
#SBATCH --error=log/interpro_%x_%j.err

set -euo pipefail

PROJECT_ROOT="/path/to/your/project"   # set this to your project root
WORKING_DIR="${PROJECT_ROOT}/GO_Terms/interproscan"
DATA_DIR="${WORKING_DIR}/data/interpro_data"
LOG_DIR="${WORKING_DIR}/log"

SHORT_FASTA="${WORKING_DIR}/data/proteomes/short_clean.faa"
LONG_FASTA="${WORKING_DIR}/data/proteomes/long_clean.faa"

mkdir -p "${LOG_DIR}" "${WORKING_DIR}/results" "${WORKING_DIR}/work"

module load nextflow/25.04.6
module load apptainer/1.4.0-rc.1

export NXF_APPTAINER_CACHEDIR="${PROJECT_ROOT}/.apptainer_cache"
export NXF_HOME="${PROJECT_ROOT}/.nextflow"

cd "${WORKING_DIR}"

echo "Running InterProScan short proteome"
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \
  -profile slurm,apptainer \
  --input "${SHORT_FASTA}" \
  --datadir "${DATA_DIR}" \
  --outdir "${WORKING_DIR}/results/medfly_short" \
  -work-dir "${WORKING_DIR}/work/short" \
  --formats tsv \
  --goterms \
  -resume

echo "Running InterProScan long proteome"
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \
  -profile slurm,apptainer \
  --input "${LONG_FASTA}" \
  --datadir "${DATA_DIR}" \
  --outdir "${WORKING_DIR}/results/medfly_long" \
  -work-dir "${WORKING_DIR}/work/long" \
  --formats tsv \
  --goterms \
  -resume

echo "All InterProScan runs completed successfully."
```

```bash
sbatch ${PROJECT_ROOT}/GO_Terms/interproscan/scripts/run_interproscan_long_short.sh
```

## 6. Convert InterProScan TSV output to a gene2go table

InterProScan's TSV output isn't in a convenient `gene_id -> GO_id` shape, so
a small Python script (`interproscan_tsv_to_gene2go.py`) extracts every
`GO:XXXXXXX` accession from each row (searching all columns, since the
column layout can vary between InterProScan member-database outputs) and
writes a two-column `source_id \t go_id` table:

```python
#!/usr/bin/env python3
import re
import sys
import csv
from collections import defaultdict

GO_RE = re.compile(r"GO:\d+")

def parse_interpro_tsv(tsv_file, out_file):
    gene2go = defaultdict(set)

    with open(tsv_file, "r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue

            seq_id = row[0].strip()
            if not seq_id:
                continue

            gos = set()
            for field in row:
                if not field:
                    continue
                gos.update(GO_RE.findall(field))

            for go_id in gos:
                gene2go[seq_id].add(go_id)

    with open(out_file, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["source_id", "go_id"])
        for seq_id in sorted(gene2go):
            for go_id in sorted(gene2go[seq_id]):
                writer.writerow([seq_id, go_id])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage: python interproscan_tsv_to_gene2go.py input.tsv output.tsv\n"
        )
        sys.exit(1)

    parse_interpro_tsv(sys.argv[1], sys.argv[2])
```

Run it on each proteome's InterProScan output:

```bash
python interproscan_tsv_to_gene2go.py \
  ${PROJECT_ROOT}/GO_Terms/interproscan/results/medfly_short/short_clean.faa.tsv \
  ${PROJECT_ROOT}/GO_Terms/interproscan/results/medfly_short/medfly_short_gene2go.tsv

python interproscan_tsv_to_gene2go.py \
  ${PROJECT_ROOT}/GO_Terms/interproscan/results/medfly_long/long_clean.faa.tsv \
  ${PROJECT_ROOT}/GO_Terms/interproscan/results/medfly_long/medfly_long_gene2go.tsv
```
