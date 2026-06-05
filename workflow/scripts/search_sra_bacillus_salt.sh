#!/usr/bin/env bash
set -euo pipefail

mkdir -p metadata/sra_tables logs/download

QUERY='((Bacillus[Organism]) AND ("RNA-Seq"[Strategy]) AND (salt OR NaCl OR salinity OR osmotic))'

echo "Searching SRA with query:"
echo "$QUERY"

esearch -db sra -query "$QUERY" \
  | efetch -format runinfo \
  > metadata/sra_tables/bacillus_salt_rnaseq_runinfo.csv

echo "Saved:"
echo "metadata/sra_tables/bacillus_salt_rnaseq_runinfo.csv"
