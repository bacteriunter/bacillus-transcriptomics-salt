#!/usr/bin/env bash
set -euo pipefail

mkdir -p raw_data/fastq logs/download

while read -r srr
do
    echo "Downloading $srr"

    fasterq-dump "$srr" \
        --split-files \
        --threads 2 \
        --outdir raw_data/fastq \
        > logs/download/${srr}.log 2>&1

    gzip -f raw_data/fastq/${srr}_*.fastq

done < metadata/download_lists/srr_accessions.txt
