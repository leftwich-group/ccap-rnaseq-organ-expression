#!/bin/bash

PROJECT_ROOT="/path/to/your/project"   # set this to your project root
RAW_DIR="${PROJECT_ROOT}/Bulk_RNA/raw_data/D:X204SC22083995-Z01-F001/01.RawData"
OUT="samplesheet.csv"
STRAND="unstranded"

echo "sample,fastq_1,fastq_2,strandedness" > $OUT

for dir in ${RAW_DIR}/RN*/; do

    sample=$(basename $dir)

    r1=$(realpath ${dir}${sample}_1.fq.gz)
    r2=$(realpath ${dir}${sample}_2.fq.gz)

    echo "${sample},${r1},${r2},${STRAND}" >> $OUT

done

echo "Samplesheet created:"
cat $OUT