# Interactive RNA-seq walkthrough (one sample)

**Author:** Jean-Charles de Coriolis — 11/3/26

**Purpose:** Look at the major conceptual steps of an RNA-seq pipeline by
running them manually on a single paired-end sample. This mirrors what
nf-core/rnaseq does internally but exposes each step.

**Major steps:**

1. Inspect raw sequencing reads
2. Quality control (FastQC)
3. Adapter / quality trimming (fastp)
4. Genome indexing (STAR)
5. Read alignment
6. Alignment statistics
7. Gene-level counting

**Final result:** A gene count table used for differential expression
analysis.

---

## 1. Log on to HPC

```bash
# on a Mac
ssh abc12def@hali.uea.ac.uk
# ~ enter password ~

# on windows
# open PuTTY
# domain : @hali.uea.ac.uk
# port   : 22
# open session
# enter UEA username (i.e. abc12def)
# enter password
```

## 2. Start an interactive HPC session

We request compute resources from the SLURM scheduler. This gives us CPUs
and RAM to run commands interactively.

```bash
srun --partition=compute --cpus-per-task=16 --mem=64G --time=04:00:00 --pty bash
```

## 3. Load required software

Bioinformatics tools are typically installed via modules or conda. We load
the tools needed for this workflow.

```bash
module purge

module load fastqc
module load fastp
module load STAR
module load samtools
module load python/anaconda/2023.07/3.11.

# Need to create a conda environment to download the subread package
conda create -n rnaseq -c bioconda subread
source activate rnaseq

# Confirm that tools are available
which fastqc
which STAR
which samtools
which featureCounts
```

## 4. Create a project directory and download data

In bioinformatics projects it is important to organise files clearly. We
will create a project directory containing:

```
project/
├── data        (RNA-seq reads)
├── genome      (reference genome + annotation)
├── results     (analysis outputs)
```

```bash
mkdir -p rnaseq_project
cd rnaseq_project

mkdir data
mkdir genome
mkdir results
```

### Download example RNA-seq data

Here we use a single paired-end sample hosted on GitHub. Replace the URL
with the repository you plan to use.

```bash
cd data

wget https://raw.githubusercontent.com/YOUR_REPOSITORY/RN1_1.fq.gz
wget https://raw.githubusercontent.com/YOUR_REPOSITORY/RN1_2.fq.gz

cd ..
```

### Download the reference genome

We download the *Ceratitis capitata* reference genome (Ccap 2.1) from the
NCBI genome FTP.

```bash
cd genome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/347/755/GCF_000347755.3_Ccap_2.1/GCF_000347755.3_Ccap_2.1_genomic.fna.gz
```

### Download the gene annotation

Gene annotations describe where genes and exons are located in the genome.
They are required for assigning reads to genes later.

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/347/755/GCF_000347755.3_Ccap_2.1/GCF_000347755.3_Ccap_2.1_genomic.gff.gz
```

### Convert GFF → GTF (required by many RNA-seq tools)

Some RNA-seq tools require the annotation in GTF format. We convert the
downloaded GFF file.

```bash
module load gffread

gffread \
GCF_000347755.3_Ccap_2.1_genomic.gff.gz \
-T \
-o genome_annotation.gtf
```

## 5. Define file path variables

To avoid writing long file paths repeatedly, we define variables.

```bash
PROJECT=$(pwd)

R1=$PROJECT/data/RN1_1.fq.gz
R2=$PROJECT/data/RN1_2.fq.gz

GENOME=$PROJECT/genome/GCF_000347755.3_Ccap_2.1_genomic.fna.gz
GTF=$PROJECT/genome/genome_annotation.gtf
```

### Verify files exist

```bash
ls data
ls genome
```

## 6. Inspect raw reads

FASTQ files contain sequencing reads. Each read consists of 4 lines:

1. read identifier
2. nucleotide sequence
3. separator (`+`)
4. quality scores

```bash
zcat $R1 | head
```

Example structure:

```
@READ_ID
ACGTTGATCGATCGATCGATC
+
FFFFFFFFFFFFFFFFFFFFF
```

Quality scores encode sequencing confidence.

## 7. Raw read quality control (FastQC)

FastQC summarises sequencing quality across the dataset.

```bash
mkdir fastqc_raw

fastqc \
--threads 16 \
--outdir fastqc_raw \
$R1 $R2
```

Output files:

```bash
ls fastqc_raw
```

Open the HTML report and examine:

- Per base sequence quality
- Adapter contamination
- GC distribution
- Overrepresented sequences

## 8. Trim reads (fastp)

Trimming removes sequencing adapters and low-quality bases.

```bash
mkdir trimmed_reads

fastp \
--in1 $R1 \
--in2 $R2 \
--out1 trimmed_reads/RN1_R1.trimmed.fq.gz \
--out2 trimmed_reads/RN1_R2.trimmed.fq.gz \
--thread 16 \
--html trimmed_reads/fastp_report.html \
--json trimmed_reads/fastp_report.json
```

Outputs:

- trimmed FASTQ files
- HTML report describing trimming statistics

## 9. Build genome index (STAR)

Aligners cannot efficiently map reads directly to a raw genome. They build a
searchable index first.

```bash
mkdir star_index

STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir star_index \
--genomeFastaFiles $GENOME \
--sjdbGTFfile $GTF \
--sjdbOverhang 149
```

Explanation: `sjdbOverhang` = read length − 1. Required for accurate splice
junction mapping.

## 10. Align reads to the genome

RNA-seq reads frequently span exon-exon junctions. STAR performs
splice-aware alignment.

```bash
mkdir star_alignment

STAR \
--runThreadN 16 \
--genomeDir star_index \
--readFilesIn trimmed_reads/RN1_R1.trimmed.fq.gz trimmed_reads/RN1_R2.trimmed.fq.gz \
--readFilesCommand zcat \
--sjdbGTFfile $GTF \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix star_alignment/RN1.
```

Main output: `RN1.Aligned.sortedByCoord.out.bam`

BAM = binary alignment file containing mapped reads.

## 11. Inspect alignment statistics

STAR produces a detailed alignment summary.

```bash
cat star_alignment/RN1.Log.final.out
```

Important metrics to discuss:

- Uniquely mapped reads
- Multi-mapped reads
- Unmapped reads
- Splice junction counts

Typical RNA-seq mapping rates: 70–90% uniquely mapped reads.

## 12. Examine the BAM file

Count total alignments:

```bash
samtools view -c star_alignment/RN1.Aligned.sortedByCoord.out.bam
```

Create BAM index (required for many downstream tools):

```bash
samtools index star_alignment/RN1.Aligned.sortedByCoord.out.bam
```

Alignment summary statistics:

```bash
samtools flagstat star_alignment/RN1.Aligned.sortedByCoord.out.bam
```

## 13. Count reads per gene

Convert aligned reads into gene counts.

```bash
mkdir counts

featureCounts \
-T 16 \
-p \
-t exon \
-g gene_id \
-a $GTF \
-o counts/RN1_counts.txt \
star_alignment/RN1.Aligned.sortedByCoord.out.bam
```

Explanation:

- `-p` paired-end data
- `-t exon` count reads overlapping exons
- `-g gene_id` summarise counts per gene

## 14. Inspect the count table

```bash
head counts/RN1_counts.txt
```

Structure:

```
GeneID    Length   RN1
gene1     1200     34
gene2     500      2
```

Each number represents the number of fragments assigned to that gene.

## 15. Connection to nf-core/rnaseq

nf-core/rnaseq automates this entire workflow and adds many additional
steps, including:

- strandedness detection
- rRNA contamination estimation
- insert size distribution
- transcript quantification
- MultiQC reporting
- automated job parallelisation

Conceptually, however, the workflow remains:

```
FASTQ
 ↓
QC
 ↓
TRIM
 ↓
ALIGN
 ↓
BAM QC
 ↓
GENE COUNTS
 ↓
DIFFERENTIAL EXPRESSION ANALYSIS
```
